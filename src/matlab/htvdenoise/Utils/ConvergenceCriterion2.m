function [OK,threshold]=ConvergenceCriterion2(f,criterion,PixelNum,record_max,branch)
% PixelNum=2048*2048;
% criterion=0.7;
if branch==1
    f(f<0)=0;
    n_f = zeros(size(f), 'single');
    n_f = ( f - min( f(:) ) ) / ( max(f(:)) - min(f(:)) );
    final_f = uint16(record_max*(n_f));
    
    threshold1 = sum(sum(final_f==0))/PixelNum;
    threshold2 = sum(sum(final_f==1))/PixelNum;
    threshold3 = sum(sum(final_f==2))/PixelNum;

    threshold = threshold1+ threshold2 + threshold3;
    OK = (threshold>=criterion);
elseif branch==0
    threshold1 = sum(sum(f==0))/PixelNum;
    threshold2 = sum(sum(f==1))/PixelNum;
    threshold3 = sum(sum(f==2))/PixelNum;

    threshold = threshold1+ threshold2 + threshold3;
    OK = (threshold>=criterion);
end

end

