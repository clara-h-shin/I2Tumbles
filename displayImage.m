function [im1,objsthisframe]=displayImage(A, im1, ...
                        showthreshprocess, showprocess, ...
                        threshprocessA, processedA, displayRange, bitDepth, ...
                        trackthisdone, islinkdone, ...
                        currframe, frmin, objs_link, objs, ...
                        hdispcircles, hdispIDs, hdisptracks, ...
                        orientstr, sizestr)

% Display the present frame
if exist('im1', 'var')
    delete(im1); clear im1
end

if showthreshprocess
    im1 = imshow(threshprocessA, []);

elseif showprocess
    im1 = imshow(processedA, []);

else
    im1 = imshow(A, []);

end

if ~isempty(displayRange)
    % Adjust display range
    caxis(displayRange*(2^bitDepth-1));
end
hold on

% display track information on the image, if desired
if trackthisdone(currframe-frmin+1)
    % tracking has been done
    if islinkdone
        objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
    else
        objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
    end
    if hdispcircles
        % plot circles on top, if tracking is done
        plot(objsthisframe(1,:), objsthisframe(2,:), 'o', 'color', [0.3 0.7 0.5])
        % If we're determining orientation, draw a line on each
        % object of length = major axis length and orientation =
        % found orientation
        if ~strcmp(orientstr, 'none')
            % make an array of line starting and ending points
            startx = objsthisframe(1,:) - objsthisframe(8,:).*cos(objsthisframe(7,:));
            starty = objsthisframe(2,:) - objsthisframe(8,:).*sin(objsthisframe(7,:));
            endx = objsthisframe(1,:) + objsthisframe(8,:).*cos(objsthisframe(7,:));
            endy = objsthisframe(2,:) + objsthisframe(8,:).*sin(objsthisframe(7,:));
            % draw lines
            plot([startx; endx], [starty; endy], '-', 'color', [0.6 1.0 0.3])
        end
        if ~strcmp(sizestr, 'none')
            %make array of lines of length equal to the object diameter
            startx = objsthisframe(1,:) - objsthisframe(7,:);
            starty = objsthisframe(2,:);
            endx = objsthisframe(1,:) + objsthisframe(7,:);
            endy = objsthisframe(2,:);
            % draw lines
            plot([startx; endx], [starty; endy], '-', 'color', [0.6 1.0 0.3])
            for k = 1:length(objsthisframe(1,:))
                if ~isnan(objsthisframe(7,k))
                    rectangle('Position',[objsthisframe(1,k)-objsthisframe(7,k) objsthisframe(2,k)-objsthisframe(7,k) 2*objsthisframe(7,k) 2*objsthisframe(7,k)], 'Curvature', [1,1],'EdgeColor', [0.6 1.0 0.3]);
                end
            end
        end
    end
    if hdispIDs
        % show IDs
        if islinkdone
            for j=1:size(objsthisframe,2)
                text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(6,j)), 'color', [1.0 0.4 0.1])
            end
        else
            for j=1:size(objsthisframe,2)
                text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(4,j)), 'color', [1.0 0.4 0.1])
            end
        end
    end
    if hdisptracks
        % plot lines corresponding to the tracks of each object
        if islinkdone
            % necessary
            for j=1:size(objsthisframe,2)
                allx = objs_link(1,objs_link(6,:)==objsthisframe(6,j));
                ally = objs_link(2,objs_link(6,:)==objsthisframe(6,j));
                plot(allx, ally, '-', 'color', [1 0.3 0.5])
            end
        end
    end
end