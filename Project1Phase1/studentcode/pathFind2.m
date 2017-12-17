function [path] = pathFind2(input,setX,setY)

%pathFind is a function that plots a path from the origin (0,0)to either a
%random or a predetermined destination. The numbers represent the movement
%in units of centimeter. (example: the coordinate (123,45) is 123 cm to the
%right and 45 cm below the origin)


%The function takes either 1 or 2 inputs. The first input determines how
%many "stop points" you want along the path; you can either input an
%integer, or type in 'R', and the function will randomly generate between
%5-10 stop points.



    
      
    

if input=='R'
    
    isreal(setX) && (floor(setX)== setX);
    isreal(setY) && (floor(setY)== setY);
    
    n=floor(6*rand()+5);
    
    leftx = floor(400*rand(1,n));
    leftxPadded = zeros(1,n+1);
    leftxPadded(1:n)=leftx;
    leftxPadded(1,n+1)=setX;
    
    
    rightyRand = floor(setY*rand(1,n));
    rightySort = sort(rightyRand);
    rightyPadded=zeros(1,n+1);
    rightyPadded(1:n)=rightySort;
    rightyPadded(1,n+1)=setY;
    path = [leftxPadded;rightyPadded]';
   

elseif isreal(input) && (floor(input)== input)
    
    isreal(setX) && (floor(setX)== setX);
    isreal(setY) && (floor(setY)== setY);

    leftx = floor(400*rand(1,input-1));
    leftxPadded = zeros(1,input);
    leftxPadded(1:input-1)=leftx;
    leftxPadded(1,input)=setX;
    
    righty = floor(setY*rand(1,input-1));
    rightySort = sort(righty);
    rightyPadded=zeros(1,input);
    rightyPadded(1:(input-1))=rightySort;
    rightyPadded(1,input)=setY;
    path = [leftxPadded;rightyPadded]';


end

%plot(path(:,1),path(:,2))
%set(gca,'Ydir','reverse')


end

