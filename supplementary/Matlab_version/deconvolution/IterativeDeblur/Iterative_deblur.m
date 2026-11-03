function data_de=Iterative_deblur(data,kernel,iteration,rule)

if nargin < 3 || isempty(iteration)
    iteration=10;
end
if nargin < 4 || isempty(rule)
    rule=1;
end

data=data./max(data(:));
if ndims(data)==3
    for i=1:size(data,3)
        data_de(:,:,i)=real(deblur_core(data(:,:,i),kernel,iteration,rule));
    end
else
    data_de=real(deblur_core(data,kernel,iteration,rule));
end