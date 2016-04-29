function expandedConfigs = ExpandConfigsRandom(configs,steps,level,npoints,deltaFact,bounds,marker_w,marker_h)

  fact = deltaFact^level;
  halfstep_tx = steps.tx/fact;
  halfstep_ty = steps.ty/fact;
  halfstep_tz = steps.tz/fact;
  halfstep_rx = steps.rx/fact;
  halfstep_rz0 = steps.rz0/fact;
  halfstep_rz1 = steps.rz1/fact;
  
  numConfigs = size(configs,1);
  
  randvec = floor(3*rand(npoints*numConfigs,6)-1); % random vectors in {-1,0,1}^6
  expanded= repmat(configs,[npoints,1]);
  ranges = [halfstep_tx,halfstep_ty,halfstep_tz,halfstep_rx,halfstep_rz0,halfstep_rz1];
  addvec = repmat(ranges,[npoints*numConfigs,1]);
  
  addvec(:,3) = (randvec(:,3) >= 0) .* addvec(:,3).*(expanded(:,3).^2) ./ (1 - abs(addvec(:,3)).*expanded(:,3))...
              + (randvec(:,3) < 0) .* addvec(:,3).*(expanded(:,3).^2) ./ (1 + abs(addvec(:,3)).*expanded(:,3));
              
  weight = expanded(:,3) + norm([marker_w, marker_h]) .* sin(expanded(:,4));
  
  % cal delta according to nearby pose
  addvec(:,1) = addvec(:,1).*weight;
  addvec(:,2) = addvec(:,2).*weight;
  positve4  = 1./addvec(:,4).*( asin(2 - 1./ (1./(2-sin(expanded(:,4))) + addvec(:,4))) - expanded(:,4));
  negative4 = 1./addvec(:,4).*( expanded(:,4) - asin(2 - 1./ (1./(2-sin(expanded(:,4))) - addvec(:,4))));
  addvec(:,4) = addvec(:,4).*((randvec(:,4) >= 0) .* positve4 + (randvec(:,4) < 0).* negative4);
  addvec(:,5) = addvec(:,5);
  addvec(:,6) = addvec(:,6);
  
  % expand poses
  expandedConfigs = expanded + randvec.*addvec;
  
  % delete unwanted poses
  [row, ~] = find(imag(expandedConfigs));
  expandedConfigs = expandedConfigs(~ismember(1:size(expandedConfigs,1), row), :);
  [row, ~] = find(expandedConfigs(:,3) > bounds.tz(2) | expandedConfigs(:,3) < bounds.tz(1));
  expandedConfigs = expandedConfigs(~ismember(1:size(expandedConfigs,1), row), :);
  [row, ~] = find(expandedConfigs(:,4) > bounds.rx(2) | expandedConfigs(:,4) < bounds.rx(1));
  expandedConfigs = expandedConfigs(~ismember(1:size(expandedConfigs,1), row), :);

end




