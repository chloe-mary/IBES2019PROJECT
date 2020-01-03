function [ECG,pulse] = liremat(path)
    load(path);
    ECG=data(:,1);
    pulse=data(:,2);
end

