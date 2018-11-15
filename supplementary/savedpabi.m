%Func
j=[11,12,13,14,15,17,19,20,21,22,23,27,28,29,30];
j=[5, 6, 7];
for i=1:size(j,2)
    temp=sprintf('H:\\OpenfMRI-ds157\\PAPI2\\FunImg\\Subject0%d',j(i));
    mkdir (temp)
    temp1=sprintf('H:\\OpenfMRI-ds157\\ds157\\sub-%d\\func\\bold.nii.gz',j(i));
    temp2=sprintf('H:\\OpenfMRI-ds157\\PAPI2\\FunImg\\Subject0%d\\',j(i));
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(temp1,temp2);
    temp=sprintf('H:\\OpenfMRI-ds157\\PAPI2\\T1Img\\Subject0%d',j(i));
    mkdir (temp)
    temp1=sprintf('H:\\OpenfMRI-ds157\\ds157\\sub-%d\\anat\\T1w.nii.gz',j(i));
    temp2=sprintf('H:\\OpenfMRI-ds157\\PAPI2\\T1Img\\Subject0%d\\',j(i));
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(temp1,temp2);
end