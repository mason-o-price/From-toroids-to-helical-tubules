function c = videcolors(N)
%COLORPALETTE 
%generate color palette by linear interpolation

colors(1,:) = [0,129,167];       %1 -> darker blue
colors(2,:) = [0, 175, 185];    %2 -> pastel blue
colors(3,:) = [253, 252, 220];   %3 -> cream
colors(4,:) = [254, 217, 183];   %4 -> tan/light salmon
colors(5,:) = [240, 113, 103];   %5 -> salmon

%output red if there is only one color
if N == 1
    c = colors(5,:);
else
    % fractional values or each given colors
    cLens = size(colors,1);
    color_section = (0:cLens-1)/(cLens-1);
    
    % fractional values of colors we want for output
    position_section = (0:N-1)/(N-1);
    
    for i = 1 : N
        %if the position_section matches color_section, just use the
        %matching color
        if sum(color_section == position_section(i)) == 1
            color_index = color_section == position_section(i);
            c(i,:) = colors(color_index,:);
            
        else
            %otherwise, interpolate between colors
            
            %the neighboring colors
            left_index = max(find(color_section < position_section(i)));
            right_index = min(find(color_section > position_section(i)));
            
            %and the distances to these colors
            left_dist = position_section(i) - color_section(left_index);
            right_dist = color_section(right_index) - position_section(i);
            
            %ratio of mixing of the two colors , left and right
            left_ratio = right_dist/(left_dist + right_dist);
            right_ratio = left_dist/(left_dist + right_dist);
            
            %finally, the color is obtained from linear interpolation
            c(i,:) = colors(left_index,:).*left_ratio + ...
                        colors(right_index,:).*right_ratio;
        end
    end
end

c = c/255;

end

