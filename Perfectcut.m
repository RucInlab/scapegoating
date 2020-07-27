function [lm] =  Perfectcut(R,k)
    lm=[];
    path = zeros(1,size(R, 1));
    pathcount = size(R, 1);
    lnotm = zeros(1,size(R, 2));
    ltrav =  zeros(1,size(R, 2));
    flag=0;
    lmadd=0;
    lnotm(k)=1;
    counttrav = 0;
   for i=1:size(R,1)
       if R(i,k)==0
           path(i)=2;
           pathcount=pathcount-1;
       end
   end
   for i=1:size(R,1)
       if path(i)==2
            for j=1:size(R,2)
                if R(i,j)==1
                   lnotm(j)=1; 
                end
            end
       end
   end
 
while pathcount>0 
   ltrav =  zeros(1,size(R, 2));
    flag=0;
    lmadd=0;
   for j = 1:size(R,2)
       if lnotm(j)==0
            for i = 1:size(R,1)
                if path(i)==0
                    if R(i,j)==1
                        ltrav(j) = ltrav(j)+1;
                    end
               end
           end
       end
       if ltrav(j)>lmadd
           flag=j;
           lmadd=ltrav(j);
       end
   end
   if flag>0
       for i=1:size(R,1)
           if path(i)==0
                if R(i,flag)==1
                    path(i)=2;
                    pathcount=pathcount-1;
                end
           end
       end
       lm=[lm,flag];
   end
   counttrav=counttrav+1;
   if counttrav > size(R,1)*2
       fprintf('No available manipulated link set\n');
       lm=[];
       break;
   end
end 
end