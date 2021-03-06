%tilecut.m
%A function to cut a stack of images into tiles

% Specify parameters here

tic;
srcnametemplate='%04d____z%d.0.tif';
targetnametemplate='img%03d_x%d_y%d.tif';
firstimgnr=3;
lastimgnr=1287;
nroftilesx=10;
nroftilesy=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imgname=sprintf(srcnametemplate,1,0);
d=dir('/mnt/data/Ashwin/ZfishExported-01072016/');
if (size(d,1)>0)
  disp(sprintf('Reading %s...',d(3).name));
  img=imread(d(3).name);
  xsize=size(img,2);
  ysize=size(img,1);
  if (xsize==0)||(ysize==0)
    disp('Image empty.');
  else
    disp(sprintf('Image size is (%d,%d).',xsize,ysize));
  
    tilesizex=floor(xsize/nroftilesx);
    tilesizey=floor(ysize/nroftilesy);
    
    %if (tilesizex*nroftilesx~=xsize)||(tilesizey*nroftilesy~=ysize)
     % disp(sprintf('ERROR: (%d,%d) is not divisible into (%d,%d) tiles!',xsize,ysize,nroftilesx,nroftilesy));
    %else
      disp(sprintf('Cutting images into %dx%d tiles of %dx%d pixels each.',nroftilesx,nroftilesy,tilesizex,tilesizey));
      for i=firstimgnr:1:lastimgnr
        imgname=d(i).name;
        d2=dir(imgname);
        if (size(d2,1)>0)
          disp(sprintf('Reading %s (%s)...',d(i).name,imgname));
          img=imread(d2(1).name);
          if (size(img,2)~=xsize)||(size(img,1)~=ysize)
            disp(sprintf('Image %s has a different size! Skipping...',d2(1).name));
          else
            for y=1:1:nroftilesy
              for x=1:1:nroftilesx
                targetdir = '/mnt/data/Ashwin/ZFishExportedChopped';
              	targetfilename=fullfile(targetdir,sprintf(targetnametemplate,i-2,x,y));
                disp(sprintf('Writing %s...',targetfilename));
                imwrite(img(((y-1)*tilesizey+1):(y*tilesizey),((x-1)*tilesizex+1):(x*tilesizex)),targetfilename);
              end;
            end;
          end;
        else
          disp(sprintf('%s not found.',imgname));
        end;
      end;
    %end;
  end;
else
  disp(sprintf('Image %s not found.',imgname));
end;
disp(sprintf('Elapsed time is %3d',toc));
