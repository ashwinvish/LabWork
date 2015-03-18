fin = open("/mnt/data0/ashwin/07122012/TrakEMFiles-S2-W004/Test/S2-W004-AffineMontage-AffineZ-ElasticStitch-ElasticZ.xml")
fout = open("/mnt/data0/ashwin/07122012/TrakEMFiles-S2-W004/Test/S2-W004-AffineMontage-AffineZ-ElasticStitch-ElasticZ.xml.tmp", "wt")
for line in fin:
    fout.write( line.replace('file_path="/data/Ashwin/', 'file_path="/mnt/data0/ashwin/') )
fin.close()
fout.close()