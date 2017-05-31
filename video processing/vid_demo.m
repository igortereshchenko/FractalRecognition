trafficVid = VideoReader('C:\studying\1_Texture_segmentation\videos\DriveByClouds.avi');
%implay('D:\Users\crystal_ship\Downloads\test_avi.avi');

nframes = 20;%trafficVid.NumberOfFrames;
I = read(trafficVid, 1);
for k = 1 : 1
    singleFrame = read(trafficVid, k);
    
    
    imwrite (singleFrame,'C:\studying\1_Texture_segmentation\videos\test.jpg');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    


FileName = 'C:\studying\1_Texture_segmentation\videos\test.jpg';

 image_origine=imread(FileName);

    if (ndims(image_origine)==3)
        image_origine=rgb2gray(image_origine);
    end

pas = 40;
Algo=1;
DrawSquares=true;
Threshold=1.2;

    image_edge=edge(image_origine,'log');

    if( DrawSquares)
        image_squares = image_origine;
    end
    
color = 0; 
pas_step = 20;
            
   MinSize=1/pas_step;
    MaxSize=1/2;
    NbBox=6;
    ShowInfos=0;
    DejaEdge=1;
    
    
    nX = floor(size(image_edge,2)/pas_step-1);
    nY = floor(size(image_edge,1)/pas_step-1);
    
    for x=0:nX
        for y=0:nY                
            if Algo
                averslope = BoxCountingFromImage(image_edge((1+y*pas_step):(y+1)*pas_step,(1+x*pas_step):(x+1)*pas_step),MinSize,MaxSize,NbBox,ShowInfos,DejaEdge);
            else
                averslope = hausdorff(image_edge((1+y*pas_step):(y+1)*pas_step,(1+x*pas_step):(x+1)*pas_step));
            end
            averslope = abs(averslope);
            M(1+y,1+x) = averslope;
            

            
                if DrawSquares && (M(1+y,1+x) > Threshold)
                    image_squares( 1+y*pas_step , (1+x*pas_step):(x+1)*pas_step )=color;
                    image_squares( (1+y)*pas_step , (1+x*pas_step):(x+1)*pas_step )=color;
                    image_squares( (1+y*pas_step):(y+1)*pas_step , 1+x*pas_step )=color;
                    image_squares( (1+y*pas_step):(y+1)*pas_step , (x+1)*pas_step )=color; 
                end;
            
        end
    end
    
    
        
        min_m=min(min(M));
        max_m=max(max(M));
        d=255.0/(max_m-min_m);
        M2=(floor((M-min_m)*d));
 figure ('Visible','off');        
       
 imagesc(M2);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) ;% Make the axes occupy the hole figure

grid_file = '';

for i=1:size(FileName,2)-4 
    grid_file(i)= FileName(i)
end;

grid_file=strcat(grid_file,'_grid.jpg');

saveas(gcf,grid_file);



%%%
close;
%%%
 figure ('Visible','off');
  %%%%figure;
  
        [mM,nM] = size(M);
        for i=1:mM
            yp(i)=22-i;
        end
        im2 = contour(1:nM,yp,M(1:mM,1:nM)); 

        
        set(gca,'XTick',[]); % Remove the ticks in the x axis!
        set(gca,'YTick',[]); % Remove the ticks in the y axis
        set(gca,'Position',[0 0 1 1]) ;% Make the axes occupy the hole figure
        
        
        contour_file = '';
        for i=1:size(FileName,2)-4 
            contour_file(i)=FileName(i);
        end
        
        
        contour_file=strcat(contour_file,'_contour.jpg');

        saveas(gcf,contour_file);
%%%%%%%%
      close;  
        %%%%%%%%%%%%%1    
        grid_im = imread(grid_file);
        origine_im = imread (FileName);
        contour_im = imread(contour_file);
        
        n1 = size (grid_im,1);
        m1 = size (grid_im,2);
        n2 = size (origine_im,1);
        m2 = size (origine_im,2);
        n3 = size (contour_im,1);
        m3 = size (contour_im,2);
        
        n = min(n1,min(n2,n3));
        m = min (m1, min(m2,m3));
        
        grid = imresize(grid_im, [n m]);
        origine = imresize(origine_im, [n m]);
        cont = imresize(contour_im, [n m]);
        transparent=imlincomb(0.15, grid, 1, origine);
    %    figure;
    %    imshow(transparent) ;
     %   title('origine + colour grid');
        transparent2=imlincomb(0.3, cont, 1, transparent);
   %     figure;
   %     imshow(transparent2) ;
    %    title('origine + colour grid + contour');
   %     figure;
        transparent3=imlincomb(0.5, cont, 0.5, origine);
   %     imshow(transparent3) ;
    %    title('origine + contour');

    


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taggedCars = zeros([n m 3 nframes], class(I));
    taggedCars(:,:,:,k) = transparent3;

end;
for k = 2 : nframes
    singleFrame = read(trafficVid, k);
    
    
    imwrite (singleFrame,'D:\Users\crystal_ship\Downloads\test.jpg');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    


FileName = 'D:\Users\crystal_ship\Downloads\test.jpg';
 image_origine=imread(FileName);

    if (ndims(image_origine)==3)
        image_origine=rgb2gray(image_origine);
    end

pas = 40;
Algo=1;
DrawSquares=true;
Threshold=1.2;

    image_edge=edge(image_origine,'log');

    if( DrawSquares)
        image_squares = image_origine;
    end
    
color = 0; 
pas_step = 20;
            

    
    MinSize=1/pas_step;
    MaxSize=1/2;
    NbBox=6;
    ShowInfos=0;
    DejaEdge=1;
    
    
    nX = floor(size(image_edge,2)/pas_step-1);
    nY = floor(size(image_edge,1)/pas_step-1);
    
    for x=0:nX
        for y=0:nY                
            if Algo
                averslope = BoxCountingFromImage(image_edge((1+y*pas_step):(y+1)*pas_step,(1+x*pas_step):(x+1)*pas_step),MinSize,MaxSize,NbBox,ShowInfos,DejaEdge);
            else
                averslope = hausdorff(image_edge((1+y*pas_step):(y+1)*pas_step,(1+x*pas_step):(x+1)*pas_step));
            end
            averslope = abs(averslope);
            M(1+y,1+x) = averslope;
            

            
                if DrawSquares && (M(1+y,1+x) > Threshold)
                    image_squares( 1+y*pas_step , (1+x*pas_step):(x+1)*pas_step )=color;
                    image_squares( (1+y)*pas_step , (1+x*pas_step):(x+1)*pas_step )=color;
                    image_squares( (1+y*pas_step):(y+1)*pas_step , 1+x*pas_step )=color;
                    image_squares( (1+y*pas_step):(y+1)*pas_step , (x+1)*pas_step )=color; 
                end;
            
        end
    end
    

    

        
        min_m=min(min(M));
        max_m=max(max(M));
        d=255.0/(max_m-min_m);
        M2=(floor((M-min_m)*d));
        
 figure ('Visible','off');      
 imagesc(M2);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) ;% Make the axes occupy the hole figure

grid_file = '';

for i=1:size(FileName,2)-4 
    grid_file(i)= FileName(i)
end;


grid_file=strcat(grid_file,'_grid.jpg');

saveas(gcf,grid_file);

  %%%%%%%%%%%%%%
  close;
  %%%%%%%%%%%
        
       figure ('Visible','off');
        [mM,nM] = size(M);
        for i=1:mM
            yp(i)=22-i;
        end
        im2 = contour(1:nM,yp,M(1:mM,1:nM)); 
        
        set(gca,'XTick',[]); % Remove the ticks in the x axis!
        set(gca,'YTick',[]); % Remove the ticks in the y axis
        set(gca,'Position',[0 0 1 1]) ;% Make the axes occupy the hole figure
        
        
        contour_file = '';
        for i=1:size(FileName,2)-4 
            contour_file(i)=FileName(i);
        end
        contour_file=strcat(contour_file,'_contour.jpg');

        saveas(gcf,contour_file);
%%%%%%%%%%%
close;
%%%%%%%%%%%
        
            
        grid_im = imread(grid_file);
        origine_im = imread (FileName);
        contour_im = imread(contour_file);
        
        n1 = size (grid_im,1);
        m1 = size (grid_im,2);
        n2 = size (origine_im,1);
        m2 = size (origine_im,2);
        n3 = size (contour_im,1);
        m3 = size (contour_im,2);
        
        n = min(n1,min(n2,n3));
        m = min (m1, min(m2,m3));
        
        grid = imresize(grid_im, [n m]);
        origine = imresize(origine_im, [n m]);
        cont = imresize(contour_im, [n m]);
        transparent=imlincomb(0.15, grid, 1, origine);
    %    figure;
    %    imshow(transparent) ;
     %   title('origine + colour grid');
        transparent2=imlincomb(0.3, cont, 1, transparent);
   %     figure;
   %     imshow(transparent2) ;
    %    title('origine + colour grid + contour');
   %     figure;
        transparent3=imlincomb(0.5, cont, 0.5, origine);
   %     imshow(transparent3) ;
    %    title('origine + contour');

    


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    

    % Convert to grayscale to do morphological processing.
   % I = rgb2gray(singleFrame);

    taggedCars(:,:,:,k) = transparent3;

end
frameRate = trafficVid.FrameRate;
implay(taggedCars,frameRate);