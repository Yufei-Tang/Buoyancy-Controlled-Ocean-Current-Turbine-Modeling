function ADCP = ADCP_Init(filename, tStart, tStop, inflowSelection)

%load the file
load(filename);

ADCP_depth1 = median(AnDepthmm/1000)-RDIBin1Mid-((1:max(SerBins))-1)*RDIBinSize; %depth range
ADCP_S1 = SerMagmmpersec'/1000; %speed in m/s
ADCP_N1 = SerNmmpersec'/1000; %North Velocity in m/s
ADCP_E1 = SerEmmpersec'/1000; %East Velocity in m/s
ADCP_D1 = -SerVmmpersec'/1000; %Down Velocity in m/s *********************

%the following code removes "bad" ADCP data from the set and replaces them with NaNs
ADCP_S_filt1 = nan*ADCP_S1;
ADCP_D_filt1 = nan*ADCP_D1;
ADCP_N_filt1 = nan*ADCP_N1;
ADCP_E_filt1 = nan*ADCP_E1;
cor64=((SerC1cnt>64)+(SerC2cnt>64)+(SerC3cnt>64)+(SerC4cnt>64))';

for ct = 1:length(ADCP_depth1)
    if ct<30
        I=find(ADCP_S1(ct,:)>0 ...
        &AnDepthmm'/1000>300 ...
        &ADCP_depth1(ct)>20 ...
        &cor64(ct,:)>2 ...
        &(SerPG4(:,ct))'>50); %18m depth is somewhat arbitrary
    else
        I=find(ADCP_S1(ct,:)>0 ...
        &AnDepthmm'/1000>300 ...
        &ADCP_depth1(ct)>20 ...
        &cor64(ct,:)>2 ...
        &(SerPG4(:,ct))'>50 ...
        &max((SerEAAcnt(:,ct-4:ct)-SerEAAcnt(:,ct-6:ct-2))')<30 ...
        &min((ADCP_S1(ct-5:ct,:)-ADCP_S1(ct-6:ct-1,:)))>-0.3); %18m depth is somewhat arbitrary
    end
    ADCP_S_filt1(ct,I) = ADCP_S1(ct,I);
    ADCP_D_filt1(ct,I) = ADCP_D1(ct,I);
    ADCP_N_filt1(ct,I) = ADCP_N1(ct,I);
    ADCP_E_filt1(ct,I) = ADCP_E1(ct,I);
end

%replace any NaNs to the previous rows values (sets the top of the water
%column to the closest good value)
ADCP_N = fillmissing(ADCP_N_filt1,'previous');
ADCP_E = fillmissing(ADCP_E_filt1,'previous');
ADCP_D = fillmissing(ADCP_D_filt1,'previous');
ADCP_S = fillmissing(ADCP_S_filt1,'previous');

TF = isnan(ADCP_N);
ADCP_N(TF) = 0;

%the date range in decimal years AD
ADCP_num_date1 = datenum([SerYear,SerMon,SerDay,SerHour+5,SerMin,0*SerMin]);

%set the variables to the length of the simulation
[~,index_start] = min(abs(ADCP_num_date1'-tStart));
[~,index_stop] = min(abs(ADCP_num_date1'-tStop));

%%this section passes variables to the main code
%original data
  ADCP.N1 = ADCP_N1(:,index_start-1:index_stop);
  ADCP.E1 = ADCP_E1(:,index_start-1:index_stop);
  ADCP.D1 = ADCP_D1(:,index_start-1:index_stop);
  ADCP.S1 = ADCP_S1(:,index_start-1:index_stop);
%filtered data without NaNs and ready for interpolation
  ADCP.u = ADCP_N(:,index_start-1:index_stop);
  ADCP.v = ADCP_E(:,index_start-1:index_stop);
  ADCP.w = ADCP_D(:,index_start-1:index_stop);
  ADCP.S = ADCP_S(:,index_start-1:index_stop);
%additional variables of interest (depth, time, index values)
  ADCP.z = ADCP_depth1;
  ADCP.t = ADCP_num_date1(index_start-1:index_stop);
  ADCP.tStart = tStart;
  ADCP.tStop = tStop;
  ADCP.iStart = index_start;
  ADCP.iStop = index_stop;
  ADCP.inflowSelection = inflowSelection;