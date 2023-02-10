function [med_bound, intQR] = quartiles(values)

prc = prctile(values,[25 50 75]);
med_bound = [prc(1)-prc(2) prc(2) prc(3)-prc(2)];
intQR = iqr(values);

disp(strcat("lower: ", num2str(med_bound(1)), "     (-",num2str(100*abs(med_bound(1)/med_bound(2))),"%)"))
disp(strcat("median: ",num2str(med_bound(2))))
disp(strcat("upper: ", num2str(med_bound(1)), "     (+",num2str(100*abs(med_bound(3)/med_bound(2))),"%)"))

disp(strcat("IQR: ",num2str(intQR)))
end