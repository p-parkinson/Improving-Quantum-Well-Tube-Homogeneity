function [med_bound, intQR] = quartiles(values)
prc = prctile(values,[25 50 75]);
med_bound = [prc(1)-prc(2) prc(2) prc(3)-prc(2)];
intQR = iqr(values);
end