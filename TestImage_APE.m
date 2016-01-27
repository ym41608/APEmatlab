function ex_mat = TestImage_APE(Marker, img, in_mat, needcompile, photometricInvariance,delta,minTz,maxTz,verbose)
    
    % adding 2 subdirectories to Matlab PATH
	AddPaths
	
	if (~exist('needcompile','var'))
		needcompile = 0;
	end
    
	if (needcompile == 1)
	    % compiling the relevant Mex-files
		CompileMex;
	end
	
	[~,ex_mat,delta,~] = APE(Marker,img,in_mat,photometricInvariance,delta,minTz,maxTz,verbose);
	

    
end
