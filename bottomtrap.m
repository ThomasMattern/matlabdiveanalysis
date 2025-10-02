function [best_tri,best_tra]=bottomtrap(d)
%function [best_tri,best_tra]=lge_test(d)

% paranoia: d is column vector
% with at least two values
% d=d(:);
N=length(d);

integral_d=sum(d);
best_traI=[0 0 0]; % inner trapezoid
best_traO=[0 0 999]; % outer trapezoid, if any
best_tri=[0 0];
for i=2:N-1;
    for j=(i+1):N-1;       
        % the trick is to find the "best" norm measure ...
    
%%%%%%%%%%%%%%%%
%%% use this if you want to take into account concav parts
%
%        % area (integral) of the trapezoid
%        % (i-1)/2*d(i) + (j-i)*(d(i)+d(j))/2 + (N-j)/2*d(j);
%        tmp = 0.5*((i-1)*d(i) + (j-i)*(d(i)+d(j)) + (N-j)*d(j));
%        tmp = tmp/integral_d;
%%%%%%%%%%%%%%%%
        % only use the convex part of the data curve
        % the overlapping partial area
        t = [d(i)/(i-1)*[0:i-1], ...
             d(i)+(d(j)-d(i))/(j-i)*[1:(j-i)], ...
             d(j)-d(j)/(N-j)*[1:(N-j)] ]';
        t = d - t;
        t(t<0)=0;
        tmp = 1-sum(t)/integral_d;
%%%%%%%%%%%%%%%%
        
        if tmp <= 1
            if tmp >= best_traI(3); best_traI = [i j tmp]; end
        else
            if tmp <= best_traO(3); best_traO = [i j tmp]; end
        end
    end;
    % area of triangle
    % should also distinguish between all/convex only area, but ...
    tmp = 0.5*(N-1)*d(i);
    tmp = tmp/integral_d;
    if tmp > best_tri(2); best_tri = [i tmp]; end  
end

if diff(best_traI(1:2)) > diff(best_traO(1:2));
    best_tra=best_traI(:);
else
    best_tra=best_traO(:);
end
best_tri=best_tri(:);
