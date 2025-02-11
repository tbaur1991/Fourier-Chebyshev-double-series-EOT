function iou = iouFCDS(X,pos_ref,or_ref,a_ref,b_ref,h_ref,nu,nth,nIOUS)
    % suppress warnings
    warning('off')

    % translate and rotate reference points and get reference poly
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = [-a_ref a_ref a_ref -a_ref; -b_ref -b_ref b_ref b_ref];
    pts_ref = R*pts_ref + pos_ref(1:2);

    % get min max of shape union
    h = c1(X(5),0,'lower');
    minh_union = min(X(3)-h/2,pos_ref(3));
    maxh_union = max(X(3)+h/2,pos_ref(3)+h_ref);

    % calculate multiple IoUs
    ious = zeros(1,nIOUS);
    hs = linspace(minh_union,maxh_union,nIOUS);
    for i = 1:nIOUS
        % get shape points of estimate
        if hs(i) < X(3)-h/2 || hs(i) > X(3)+h/2
            pts = [0;0];
        else
            hh = hs(i) - X(3);
            u = 2*hh/h;
            pts = getShapePoints2dTop(X,nu,nth,u);
        end

        % get shape points of reference
        if hs(i) < pos_ref(3) || hs(i) > pos_ref(3)+h_ref
            pts_r = [0;0];
        else
            pts_r = pts_ref;
        end
    
        % calculate intersection over union
        ious(i) = calculateIoU(pts, pts_r);
    end

    % mean value
    iou = mean(ious,2);
end

function pts = getShapePoints2dTop(X,nu,nth,u)
    % get pose and shape parameters
    pos = X(1:2); or = X(4); shape = X(6:end);

    % calculate shape points
    thetas = 0:2*pi/30:2*pi-2*pi/30;
    rs = zeros(1,length(thetas));
    for i = 1:length(thetas)
        rs(i) = fourier_chebychev_series(shape,thetas(i),u,nu,nth);
    end
    pts = [rs.*cos(thetas); rs.*sin(thetas)];

    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts, pts_ref)
    % polyshapes
    poly = polyshape(pts(1,:),pts(2,:));
    poly_ref = polyshape(pts_ref(1,:),pts_ref(2,:));

    % intersection area
    poly_int = intersect(poly,poly_ref);
    area_int = area(poly_int);

    % union area
    poly_union = union(poly,poly_ref);
    area_union = area(poly_union);

    % intersection over union
    iou = area_int/area_union;
end