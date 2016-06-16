
function color = colorplot(lines)
%Generates RGB values for lines in progressive color scale
%red  [1 0 0]    green [0 1 0]    blue [0 0 1]

color=zeros(lines,3);
a=floor(lines/2);
%red to green
int=1/(a-1);
color(1:a,1)=1:-int:0;
color(1:a,2)=0:int:1;


%green to blue
int2=1/(length(a+1:lines)-1);
color(a+1:lines,2)=1:-int2:0;
color(a+1:lines,3)=0:int2:1;


