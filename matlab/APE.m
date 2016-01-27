function [bestConfig,ex_mat,delta,sampledError] = APE(marker,img,in_mat,... % mandatory
   photometricInvariance,delta,minTz,maxTz,verbose,epsilon) %optional                                                                        
	
	% set default values for optional variables
	if (~exist('photometricInvariance','var'))
		photometricInvariance = 0;
	end
	if (~exist('delta','var'))
		delta = 0.25;
	end
	if (~exist('minTz','var'))
		minTz = 3;
	end
	if (~exist('maxTz','var'))
		maxTz = 8;
	end
	if (~exist('verbose','var'))
		verbose = 0;
	end
	if (~exist('epsilon','var'))
		epsilon = 0.15;
	end
	if (~exist('in_mat', 'var'))  %% get intrinsic based on camera image
		[h2,w2,~] = size(img);
		focal_length = norm([h2, w2]);
		in_mat = [focal_length,0,h2/2,0;0,-focal_length,w2/2,0;0,0,1,0;0,0,0,1];
	end
	
	% ensure the data type of images is double
	if ( ~strcmp(class(img),'double') || ~strcmp(class(marker),'double')) %#ok<STISA>
		error('img and marker should both be of class ''double'' (in the range [0,1])');
	end
	
	% preCalculation
	[marker, img, bounds, steps, dim] = preCal(in_mat, marker, img, delta, minTz, maxTz);
	
	% coarse-to-fine estimation
	[bestConfig,ex_mat,delta,sampledError] = C2Festimate(marker,img,in_mat,bounds,steps,dim,epsilon,delta,photometricInvariance,verbose);
	
	return
end
