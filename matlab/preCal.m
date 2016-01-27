function [marker, img, bounds, steps, dim] = preCal(in_mat, marker, img, delta, minTz, maxTz)
	
	% dimensions of images
	[hm, wm, ~] = size(marker);
	[hi, wi, ~] = size(img);
	dim.marker.w = wm;
	dim.marker.h = hm;
	dim.img.w = wi;
	dim.img.h = hi;
	
	% intrinsic parameter
	Sxf = in_mat(1,1);
	Syf = -in_mat(2,2);
	imgX0 = in_mat(1,3);
	imgY0 = in_mat(2,3);
	
	% search range in pose domain
	marker_w = wm/min(wm, hm)*0.5;
	marker_h = hm/min(wm, hm)*0.5;
	Tx_width = (imgX0 * maxTz) / Sxf - min(marker_w,marker_h);
	Ty_width = (imgY0 * maxTz) / Syf - min(marker_w,marker_h);
	minTx = -Tx_width;
	maxTx = Tx_width;
	minTy = -Ty_width;
	maxTy = Ty_width;
	minRx = 0;
	maxRx = 80*pi/180;
	minRz = -pi;
	maxRz = pi;
	dim.marker_w = marker_w;
	dim.marker_h = marker_h;
	
	% cal steps and bounds
	bounds.tx = [minTx, maxTx];
	bounds.ty = [minTy, maxTy];
	bounds.tz = [minTz,maxTz];
	bounds.rx = [minRx,maxRx];
	bounds.rz = [minRz,maxRz];
	mdian_tz = sqrt(minTz*maxTz);
	steps.tx = delta/sqrt(2)/mdian_tz;
	steps.ty = delta/sqrt(2)/mdian_tz;
	steps.tz = delta/sqrt(2)/mdian_tz;
	steps.rx = delta/sqrt(2)/mdian_tz;
	steps.rz0 = delta*sqrt(2)/mdian_tz;
	steps.rz1 = delta*sqrt(2)/mdian_tz;
	
	% smooth images
	blur_sigma = calSigmaValue(rgb2gray(marker), Sxf, Syf, dim, minTz*maxTz);
	fprintf('blur sigma : %d\n', blur_sigma);
	blur_size = 4 * blur_sigma + 1;
	params.blur_kernel  = fspecial('gaussian', blur_size, blur_sigma);
	marker = imfilter(marker,params.blur_kernel,'symmetric');
	img = imfilter(img,params.blur_kernel,'symmetric');
	imwrite(marker, 'marker.png');
  imwrite(img, 'img.png');
	% rgb to ycbcr
	marker = rgb2ycbcrNorm(marker);
	img = rgb2ycbcrNorm(img);
end

function img = rgb2ycbcrNorm(img)
	img = rgb2ycbcr(img);
	img(:,:,1) = (img(:,:,1) - 16/255) / (235/255-16/255);
	img(:,:,2) = (img(:,:,2) - 16/255) / (240/255-16/255);
	img(:,:,3) = (img(:,:,3) - 16/255) / (240/255-16/255);
end
