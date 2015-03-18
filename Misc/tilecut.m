%tilecut.m
%A function to cut a stack of images into tiles

% Specify parameters here
srcnametemplate='%03d__*.tif';
targetnametemplate='img%03d_x%d_y%d.tif';
firstimgnr=1;
lastimgnr=100;
nroftilesx=2;
nroftilesy=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgname=sprintf(srcnametemplate,1);
d=dir(imgname);
if (size(d,1)>0)
  disp(sprintf('Reading %s (%s)...',d(1).name,imgname));
  img=imread(d(1).name);
  xsize=size(img,2);
  ysize=size(img,1);
  if (xsize==0)||(ysize==0)
    disp('Image empty.');
  else
    disp(sprintf('Image size is (%d,%d).',xsize,ysize));
  
    tilesizex=round(xsize/nroftilesx);
    tilesizey=round(ysize/nroftilesy);
    
    if (tilesizex*nroftilesx~=xsize)||(tilesizey*nroftilesy~=ysize)
      disp(sprintf('ERROR: (%d,%d) is not divisible into (%d,%d) tiles!',xsize,ysize,nroftilesx,nroftilesy));
    else
      disp(sprintf('Cutting images into %dx%d tiles of %dx%d pixels each.',nroftilesx,nroftilesy,tilesizex,tilesizey));
      for i=firstimgnr:1:lastimgnr
        imgname=sprintf(srcnametemplate,i);
        d=dir(imgname);
        if (size(d,1)>0)
          disp(sprintf('Reading %s (%s)...',d(1).name,imgname));
          img=imread(d(1).name);
          if (size(img,2)~=xsize)||(size(img,1)~=ysize)
            disp(sprintf('Image %s has a different size! Skipping...',d(1).name));
          else
            for y=1:1:nroftilesy
              for x=1:1:nroftilesx
                targetfilename=sprintf(targetnametemplate,i,x,y);
                disp(sprintf('Writing %s...',targetfilename));
                imwrite(img(((y-1)*tilesizey+1):(y*tilesizey),((x-1)*tilesizex+1):(x*tilesizex)),targetfilename);
              end;
            end;
          end;
        else
          disp(sprintf('%s not found.',imgname));
        end;
      end;
    end;
  end;
else
  disp(sprintf('Image %s not found.',imgname));
end;