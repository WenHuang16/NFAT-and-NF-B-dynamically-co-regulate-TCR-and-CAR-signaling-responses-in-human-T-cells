%% loading_data
% NC ratio are exported here
% RAW data can be found in rls after loading nucle-to-cyto-*.mat

fls_all=[];

   for pos=1:16
       cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2021\2021-01-21\z_t0\pos' num2str(pos)])
       fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2021\2021-01-21\z_t0\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
       fls_all=[fls_all;fls_pos];
   end


fls_all_r = [];
for i= 2*(1:length(fls_all)/2)
    fls_all_r = [fls_all_r, fls_all(i)];
end

dirname = unique({fls_all_r.folder}');
fls_all_data = struct('folder',{fls_all_r.folder},'name',{fls_all_r.name},'NFKB',{[]},'NFAT',{[]},'CalRFP',{[]},'NuciRFP',{[]});


%% raw n/c ratio, CalRFP, NuciRFP
% chanel 1:iRFP; 2:NFKB, 3:NFAT; 4:RFP
for i =1:length(fls_all_r)
    cd (fls_all_r(i).folder);
    load (fls_all_r(i).name);
    %NC RATIO  2min/frame,4h(120 frames)
    size_mena_ntc = size(mean_ntc);
    fls_all_data(i).NFKB = mean_ntc(2).data (:,1:121);
    fls_all_data(i).NFAT = mean_ntc(3).data (:,1:121);
    
    fls_all_data(i).CalRFP = [];
    fls_all_data(i).NuciRFP = [];
    for f = 1:121
        fls_all_data(i).CalRFP = [fls_all_data(i).CalRFP, mean(rls(f).channel(4).nucl)];
        fls_all_data(i).NuciRFP = [fls_all_data(i).NuciRFP, mean(rls(f).channel(1).nucl)];
    end
end




%% [HAMPEL & SMOOTH] n/c ratio, CalRFP, NuciRFP
% chanel 1:iRFP; 2:NFKB, 3:NFAT; 4:RFP
processed_data = struct('folder',{fls_all_r.folder},'name',{fls_all_r.name},'NFKB',{[]},'NFAT',{[]},'CalRFP',{[]},'NuciRFP',{[]});

Ca_NFAT = [];
Ca_indicator = [];
Ca_nuc = [];

for i =1:length(fls_all_data)
        h_NFKB = hampel(fls_all_data(i).NFKB);
        order1 = 3;
        framelen1 = 11;
        processed_data(i).NFKB = sgolayfilt(h_NFKB, order1, framelen1); %NFkB

        h_NFAT = hampel(fls_all_data(i).NFAT);
        order2 = 3;
        framelen2 = 11;
        processed_data(i).NFAT = sgolayfilt(h_NFAT, order2, framelen2); %NFAT
        Ca_NFAT =[Ca_NFAT; processed_data(i).NFAT]; % export
        
        processed_data(i).CalRFP = fls_all_data(i).CalRFP; %Ca
        Ca_indicator = [Ca_indicator; processed_data(i).CalRFP]; % export

        h_NuciRFP = hampel(fls_all_data(i).NuciRFP);
        processed_data(i).NuciRFP = h_NuciRFP; %Nu     
        Ca_nuc = [Ca_nuc; processed_data(i).NuciRFP]; % export
        
        splitStr = regexp(processed_data(i).folder,'\\','split');
        processed_data(i).pos =splitStr{1, 9} ;
        
end




    

   %% Plot together wiht 2 channels
    figure;
    time = 2*(0:120);
    % 1:79
    % 80:175
    % 176:232
    for celli =176:232
        subplot(8,12,celli-175);
        xlim([0,240]);
               
        yyaxis left;
        plot (time, processed_data(celli).NFAT ,'LineWidth',2);
         ylim([0,15]);
        %ylabel('NFAT N/C ratio');
        
        yyaxis right;
        plot (time, processed_data(celli).CalRFP,'red','LineWidth',2);
        hold on;
         ylim([0,150]);
       %ylabel('CalRFP');
        
        %xlabel('Time(min)');
       title(processed_data(celli).pos);
        
    end
    


         