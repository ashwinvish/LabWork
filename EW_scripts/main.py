import httplib # catching httplib exception
import numpy as np
import requests
import json
from retrying import retry
import struct
import multiprocessing
import os
import glob
import sys
from time import sleep

def retry_if_backend_error(exception):
    if isinstance(exception, httplib.InternalServerError):
        return True
    elif isinstance(exception, httplib.IncompleteRead):
        return True
    return False

@retry(retry_on_exception=retry_if_backend_error, wait_exponential_multiplier=1000, wait_exponential_max=10000)
def download_file(url):
    response = requests.get(url)
    return response.content

def retrieve_tasks(cid):
    tasks = json.loads(download_file("http://eyewire.org/1.0/cell/{}/tasks".format(cid)))
    return tasks

def dstrip_to_ply(filename, dstrip, world_min, voxel_res):
    if (dstrip.size < 18 or dstrip.size % 6 != 0):
        raise ValueError("Incorrect dstrip length")

    dstrip = dstrip.reshape(dstrip.size / 3, 3)
    dstrip[0::2] = np.array(map(lambda v: 0.5 * v * voxel_res + world_min, dstrip[0::2]))
    dstrip = dstrip.ravel()

    vertex_count = dstrip.size / 6
    triangle_count = vertex_count - 2

    ply = open(filename, "wb")
    ply.write("ply\n")
    ply.write("format binary_little_endian 1.0\n")
    ply.write("element vertex {}\n".format(vertex_count))
    ply.write("property float x\n")
    ply.write("property float y\n")
    ply.write("property float z\n")
    ply.write("property float nx\n")
    ply.write("property float ny\n")
    ply.write("property float nz\n")
    ply.write("element face {}\n".format(triangle_count))
    ply.write("property list uchar int vertex_indices\n")
    ply.write("end_header\n")

    vertex_struct = struct.Struct('<' + 'f' * dstrip.size)
    ply.write(vertex_struct.pack(*dstrip))

    (v0, v1, v2) = (0,1,2)
    face_struct = struct.Struct('<Biii')

    for i in xrange(triangle_count):
        ply.write(face_struct.pack(3, v0, v1, v2))
        if i % 2 == 0: v0 = v2
        else:          v1 = v2
        v2 += 1

    ply.close()

def dstrip_to_obj(filename, dstrip, world_min, world_max):
    if (dstrip.size < 18 or dstrip.size % 6 != 0):
        raise ValueError("Incorrect dstrip length")

    vertex_count = dstrip.size / 6
    triangle_count = vertex_count - 2

    obj = open(filename, "w")
    for i in xrange(0, dstrip.size, 6):
        obj.write("v {} {} {}\n".format(str(dstrip[i]), str(dstrip[i+1]), str(dstrip[i+2])))
        obj.write("vn {} {} {}\n".format(str(dstrip[i+3]), str(dstrip[i+4]), str(dstrip[i+5])))

    (v0, v1, v2) = (1,2,3)

    for i in xrange(triangle_count):
        obj.write("\nf {}//{} {}//{} {}//{}".format(v0, v0, v1, v1, v2, v2))
        if i % 2 == 0: v0 = v2
        else:          v1 = v2
        v2 += 1

    obj.close()

def worker_main(queue):
    f = open("{}/tmp/{}.pid".format(sys.path[0], os.getpid()), "w")
    f.close()
    while True:
        task = queue.get(True)

        if (task is None):
            queue.put(None) # Put it back for other running workers
            os.remove("{}/tmp/{}.pid".format(sys.path[0], os.getpid()))

            print "End of queue - killing worker {}".format(os.getpid())
            return
        else:
            #cid, mesh_nr, tid, path, xmin, ymin, zmin, xmax, ymax, zmax = task

            url = "https://storage.googleapis.com/overview_meshes/meshes/{}/{}/2.dstrip".format(task["cid"], task["tid"])
            print url
            sys.stdout.flush()
            blob = download_file(url)
            data = np.fromstring(blob, dtype=np.float32)

            path = "{}/tmp/{}/{}".format(sys.path[0], task["cid"], task["mesh_nr"])

            try:
                dstrip_to_ply(path + "_dstrip.ply", data, task["world_min"], task["voxel_res"])
            except ValueError:
                print "Invalid mesh: https://storage.googleapis.com/overview_meshes/meshes/{}/{}/2.dstrip".format(task["cid"], task["tid"])
                continue

            if os.path.exists("{}.ply".format(path)):
                os.remove("{}.ply".format(path))

            cmd = "meshlabserver -i " + path + "_dstrip.ply -o " + path + ".ply -m vn -s cleanup.mlx &> /dev/null"
            os.system(cmd)

            i = 0
            while i < 300 and not os.path.exists("{}.ply".format(path)):
                sleep(0.1)
                i += 1
            if i == 300:
                print "Timeout when saving mesh {}.ply".format(path)
                continue

def main():
    if len(sys.argv) == 1:
        print "Usage main.py <cell_id_1> [<cell_id_2> ...]"
        return
    
    # Create tmp directory
    if not os.path.exists("{}/tmp".format(sys.path[0])):
        os.makedirs("{}/tmp".format(sys.path[0]))

    for arg in sys.argv:
        try:
            cid = int(arg)
        except ValueError:
            continue

        # The workers call meshlabserver, which detaches itself from the commandline, making it
        # difficult to tell whether it finished or not
        # To solve this each worker creates a dummy pid file when started, and deletes it when done.

        # Get all existing .pid files -- we only want to wait for deletion of our own ones
        oldpid = set(glob.glob('./tmp/*.pid'))

        # Create Queue and Pool
        queue = multiprocessing.Queue()
        pool = multiprocessing.Pool(8, worker_main, (queue,))
        sleep(1)

        # Store the pid files that were newly generated
        newpid = set(glob.glob('./tmp/*.pid')) - oldpid

        # Create directory for cell
        if not os.path.exists("{}/tmp/{}".format(sys.path[0], cid)):
            os.makedirs("{}/tmp/{}".format(sys.path[0], cid))

        # Get dataset and cell info
        tasks = retrieve_tasks(cid)

        # Add all tasks to the queue
        mesh_cnt = 0
        for task in tasks['tasks']:
            params = {
                "cid": cid,
                "mesh_nr": mesh_cnt,
                "tid": task['id'],
                "world_min": np.array([task['bounds']['min']['x'], task['bounds']['min']['y'], task['bounds']['min']['z']]),
                "world_max": np.array([task['bounds']['max']['x'], task['bounds']['max']['y'], task['bounds']['max']['z']]),
                "voxel_res": np.array([tasks['world']['volumes']['resolution']['x'], tasks['world']['volumes']['resolution']['y'], tasks['world']['volumes']['resolution']['z']])
            }

            queue.put(params)
            mesh_cnt += 1

        # Add poison pill to shutdown workers
        queue.put(None)

        # Wait until all workers finished and deleted their pid file (might get stuck!)
        while True:
            done = True
            for f in newpid:
                if os.path.exists("{}/tmp/{}".format(sys.path[0], f)):
                    done = False
            if done is True:
                break
            sleep(1)
        
        # Clean up
        pool.close()
        pool.join()

        # Merge all downloaded, processed meshes into one single file
        inputfiles = " ".join(map(lambda x: "-i ./tmp/{}/{}.ply".format(cid, x), range(mesh_cnt)))
        cmd = "meshlabserver " + inputfiles + " -o " + "{}/{}".format(sys.path[0], cid) + ".ply -m vn -s merge.mlx"
        os.system(cmd)

if __name__ == '__main__':
    main()

