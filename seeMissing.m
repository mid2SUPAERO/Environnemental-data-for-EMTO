parpool('local',4)
density=0.01:0.1:1;
m=size(density,2);
parfor k =1:1
    missingname='0mis';
    missingID=fopen(missingname,'w');
    for j=1:m
        locfilename=['pointHRr3',num2str(j)];
        if ~isfile(locfilename)
            fprintf(missingID,'%14.10f',j);
            fprintf(missingID,'\n');
        end
    end
end
    
fclose('all');



