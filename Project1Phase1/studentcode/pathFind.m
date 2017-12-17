function [path] = pathFind(input,varargin)

%pathFind is a function that plots a path from the origin (0,0)to either a
%random or a predetermined destination. The numbers represent the movement
%in units of centimeter. (example: the coordinate (123,45) is 123 cm to the
%right and 45 cm below the origin)

%The function has already been configured to fall within the 
%RMC arena parameters (aka robot will NOT travel outside of the arena
%space)

%The function takes either 1 or 2 inputs. The first input determines how
%many "stop points" you want along the path; you can either input an
%integer, or type in 'R', and the function will randomly generate between
%5-10 stop points.

%The second input will set the final x coordinate; the final y coordinate
%will always be 738 as this is the "end" of the arena. If setting a final x
%coordinate, make sure the integer is between 0-375. If you would like a
%random endpoint, simply leave the second input space blank. Have fun!


if isempty(varargin)
    varargin ={floor(376*rand())};
end      
      
    

if input=='R'
    
    n=floor(6*rand()+5);
    leftx = floor(376*rand(1,n));
    leftxPadded = zeros(1,n+1);
    leftxPadded(1:n)=leftx;
    leftxPadded(1,n+1)=varargin{1};
    
    rightyRand = floor(739*rand(1,n));
    rightySort = sort(rightyRand);
    rightyPadded=zeros(1,n+1);
    rightyPadded(1:n)=rightySort;
    rightyPadded(1,n+1)=738;
    path = [leftxPadded;rightyPadded]';
   

elseif isreal(input) && (floor(input)== input)

    leftx = floor(376*rand(1,input-1));
    leftxPadded = zeros(1,input);
    leftxPadded(1:input-1)=leftx;
    leftxPadded(1,input)=varargin{1};
    
    righty = floor(739*rand(1,input-1));
    rightySort = sort(righty);
    rightyPadded=zeros(1,input);
    rightyPadded(1:(input-1))=rightySort;
    rightyPadded(1,input)=738;
    path = [leftxPadded;rightyPadded]';
 

else 
    fprintf('Error: Please input R or an Integer in first input, and an integer between 0-375 in second input\n');

end

%plot(path(:,1),path(:,2))
%set(gca,'Ydir','reverse')


end

