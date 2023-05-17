%% Export all single cell data by 2nd-26# monoclone
% NC ratio are exported here
% RAW data can be found in rls after loading nucle-to-cyto-*.mat

clear all
fls_all=[];
% 
% for imagedate = [10,11,13,14,15,16,17,20,21,22]
%     for pos=1:3
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = '0.5CC';
%         end
%         fls_all=[fls_all;fls_pos];
%     end
%     
%     for pos=4:6
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = 'CC';
%         end
%         fls_all=[fls_all;fls_pos];
%     end
%     
%     for pos=7:9
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = '2CC';
%         end
%         fls_all=[fls_all;fls_pos];
% 
%     end
% end
% %%
% for imagedate = [14,15,21,22]
%     for pos=10:15
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = '4CC';
%         end
%         
%         fls_all=[fls_all;fls_pos];
% 
%     end
% end
% 
% %% PI (I1.4P50)
% for imagedate = [10,11,13,14,15,16,17,20,21,22]
%     for pos= 22:24
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = 'PI';
%         end
%         
%         fls_all=[fls_all;fls_pos];
% 
%     end
%     
%     
% end

%% P (P50), I (I1.4)
for imagedate = [10,11,13,16,17,20]
%     for pos= 10:12
%         cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
%         for kk=1:length(fls_pos)
%             fls_pos(kk).treatment = 'P';
%         end
%         
%         fls_all=[fls_all;fls_pos];
% 
%     end
   for pos= 13:15
        cd (['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos)])
        fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
        fls_pos=dir(['D:\LabNAS_Backup\HuangWen\Movie\Dragonfly\2019-09-' num2str(imagedate) '\z\pos' num2str(pos) '\nucle-to-cyto-*.mat']);
        for kk=1:length(fls_pos)
            fls_pos(kk).treatment = 'I';
        end
        
        fls_all=[fls_all;fls_pos];

    end
    
    
end

dirname = unique({fls_all.folder}');
fls_all_data = struct('folder',{fls_all.folder},'name',{fls_all.name},'treatment',{fls_all.treatment},'NFKB9h',{[]},'NFAT9h',{[]},'mCherry9h',{[]});

%%
NFKB_EXCEL = []; 
NFAT_EXCEL = []; 
mCherry_EXCEL =[]; 

NFKB_ABSNI_sum_EXCEL = [];
NFAT_ABSNI_sum_EXCEL = [];
NFKB_ABSNI_mean_EXCEL = [];
NFAT_ABSNI_mean_EXCEL = [];

for i =1:size(fls_all)
cd (fls_all(i).folder);
load (fls_all(i).name);

%Absolute Nuclear Intensity within 9h
for j= 1:108
    fls_all_data(i).NFKB_sum(j) = sum(rls(j).channel(2).nucl);
    fls_all_data(i).NFKB_mean(j) = mean(rls(j).channel(2).nucl);
    
    fls_all_data(i).NFAT_sum(j) = sum(rls(j).channel(3).nucl);
    fls_all_data(i).NFAT_mean(j) = mean(rls(j).channel(3).nucl);
end

%NC RATIO & Reporter intensity within 9h
size_mena_ntc = size(mean_ntc);
fls_all_data(i).NFKB9h = mean_ntc(2).data (:,1:108);
fls_all_data(i).NFAT9h = mean_ntc(3).data (:,1:108);
    if size_mena_ntc(2)==4
        fls_all_data(i).mCherry9h = mean_ntc(4).data (:,1:108);
    else
        fls_all_data(i).mCherry9h =linspace(0,0,108);
    end

splitStr = regexp(fls_all(i).folder,'\\','split');
fls_all_data(i).dateID = splitStr{1, 6} ;
fls_all_data(i).treatmentID =splitStr{1, 8} ;

NFKB_EXCEL = [NFKB_EXCEL;hampel(fls_all_data(i).NFKB9h)]; %%Drop outliers
NFAT_EXCEL= [NFAT_EXCEL;hampel(fls_all_data(i).NFAT9h)];
mCherry_EXCEL =[mCherry_EXCEL;hampel(fls_all_data(i).mCherry9h)];

NFKB_ABSNI_sum_EXCEL = [NFKB_ABSNI_sum_EXCEL; hampel(fls_all_data(i).NFKB_sum)];
NFKB_ABSNI_mean_EXCEL = [NFKB_ABSNI_mean_EXCEL; hampel(fls_all_data(i).NFKB_mean)];

NFAT_ABSNI_sum_EXCEL = [NFAT_ABSNI_sum_EXCEL; hampel(fls_all_data(i).NFAT_sum) ];
NFAT_ABSNI_mean_EXCEL = [NFAT_ABSNI_mean_EXCEL; hampel(fls_all_data(i).NFAT_mean)];

% scelldata(i).dataID =splitStr{1, 6};
% scelldata(i).treatmentID =splitStr{1, 8} ;

end




%% Select and smooth data

cellidx=0;
celldata=[];
time=(1:108)*5; %5min/frame,1st frame=5min


%% filter and smooth 
% no basel substrating procedure 20200521 % After substrating basel by R
NFKB_EXCEL2 = []; 
NFAT_EXCEL2 = []; 
mCherry_EXCEL2 = [];
mCherry_EXCEL3 = [];


NFKB_ABSNI_mean_EXCEL2 = []; 
NFAT_ABSNI_mean_EXCEL2 = []; 

for i =1:length(mCherry_EXCEL)    
hcelldata_NFKBi = NFKB_EXCEL(i,:);
hcelldata_NFATi = NFAT_EXCEL(i,:);
hcelldata_mCherryi = mCherry_EXCEL(i,:);

order1 = 3;
framelen1 = 11;
scelldata_NFKBi= sgolayfilt(hcelldata_NFKBi,order1,framelen1); %NFkB
NFKB_EXCEL2 = [NFKB_EXCEL2;scelldata_NFKBi] ;

order2 = 3;
framelen2 = 11;
scelldata_NFATi= sgolayfilt(hcelldata_NFATi,order2,framelen2); %NFAT
NFAT_EXCEL2 = [NFAT_EXCEL2;scelldata_NFATi] ;

scelldata_mCherryi=  (smooth(time,hcelldata_mCherryi,40,'sgolay'))'; %mCherry
mCherry_EXCEL2 = [mCherry_EXCEL2;scelldata_mCherryi] ;

Cfit = fit(time', hcelldata_mCherryi','linearinterp');
mCherry_EXCEL3 = [mCherry_EXCEL3;];


hcelldata_NFKBi_ABSNI = NFKB_ABSNI_mean_EXCEL(i,:);
hcelldata_NFATi_ABSNI = NFAT_ABSNI_mean_EXCEL(i,:);

order1 = 3;
framelen1 = 11;
scelldata_NFKBi_ABSNI= sgolayfilt(hcelldata_NFKBi_ABSNI,order1,framelen1); %NFkB_ABSNI
NFKB_ABSNI_mean_EXCEL2 = [NFKB_ABSNI_mean_EXCEL2;scelldata_NFKBi_ABSNI]; 

order1 = 3;
framelen1 = 11;
scelldata_NFATi_ABSNI= sgolayfilt(hcelldata_NFATi_ABSNI,order1,framelen1); %NFAT_ABSNI
NFAT_ABSNI_mean_EXCEL2 = [NFAT_ABSNI_mean_EXCEL2;scelldata_NFATi_ABSNI]; %NFAT_ABSNI
end


S_NFKB_ABSNI= [];
S_NFAT_ABSNI= [];
for i= 1:length(NFKB_ABSNI_mean_EXCEL)
  S_NFKB_ABSNI= [S_NFKB_ABSNI; trapz(5*(1:108),NFKB_ABSNI_mean_EXCEL2(i,:))];
  S_NFAT_ABSNI= [S_NFAT_ABSNI; trapz(5*(1:108),NFAT_ABSNI_mean_EXCEL2(i,:))];
end    
    

%%
% %% Plot together wiht 2 channels
% figure;
% for fl=1:length(fls_all)
%     subplot(6,8,fl);
%     xlim([5,540]);
%     yyaxis left;
%    
%     plot (time,scelldata(fl).sNFKB9h,'b','LineWidth',2);
%     hold on;
%     scatter( selectNFkB(fl).locs*5, selectNFkB(fl).pks,'filled','v','m','LineWidth',1);
%     ylim([0,max(scelldata(fl).sNFKB9h)*1.2]);
%    
%     % ylabel('NF-kB N/C ratio');
%     
%     yyaxis right;
%     plot (time, scelldata(fl).sNFAT9h, 'color',[50,205,50]/255,'LineWidth',2);
%     ylim([0,max(scelldata(fl).sNFAT9h)*1.2]);
% %     ylim([0,4]);
%     % ylabel('NFAT N/C ratio');
%     
%     % xlabel('Time(min)');
%     % title(['Pos' num2str(fls),'Cell'num2str(cellnumber)]);
% 
% end

% %%FOR PLOTTING EXAMPLE
% scelldata_NFKBi_ABSNI= sgolayfilt(mean_ntc(2).data ,order1,framelen1); %NFkB_ABSNI
% mean_ntc_sgolayfit = [mean_ntc_sgolayfit;scelldata_NFKBi_ABSNI]; 
% 
% scelldata_NFATi_ABSNI= sgolayfilt(mean_ntc(3).data ,order1,framelen1); %NFAT_ABSNI
% mean_ntc_sgolayfit = [mean_ntc_sgolayfit;scelldata_NFATi_ABSNI]; %NFAT_ABSNI

       
       