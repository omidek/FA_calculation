function dist = get_dist(x,y)
% calculates the distance between points x and y
% x and y sets points of m-by-2 or 3 
%   by Omid Ekrami, 2018

%     if size(x) ~= size(y)
%         msg = 'The matrices must have the same dimensions';
%         error(msg)
%     end

    dist = zeros(size(x,1),1);
    for i=1:size(x,1)
        switch size(x,2)
            case 2
                dist(i) = sqrt((x(i,1)-y(i,1))^2+(x(i,2)-y(i,2))^2);
            case 3
                dist(i) = sqrt((x(i,1)-y(i,1))^2+(x(i,2)-y(i,2))^2+(x(i,3)-y(i,3))^2);
        end
    end
end