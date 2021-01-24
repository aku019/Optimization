close all;
clear all;
clc;
%function out = PSO1_modified_1(problem, params,c)%,ncomp)


%% DECLARATION PART

nirspectra='C:\Users\MANI\Downloads\carbon-20171224T081343Z-001\carbon\all spectra\';    % Folder path where the spectroscopic samples are kept
[Refdata]=xlsread('C:\Users\MANI\Desktop\carbon-20171224T081343Z-001\carbon\Refdata_curbon.xls');%Xls sheet of polyphenol content of the 55 samples
filename_string='S%d_S%d_%d.ABS';
samp_nocompile=15;
pos_list=1;
replicant=1:5;     
spectrano_sample=(length(pos_list)*length(replicant));
startrow1=13;
endrow1=471;

All_abs=[];
%% COMPILATION OF RAW DATA PART

for sample_id=1:samp_nocompile 
for pos=1:length(pos_list)                                                    
for repl=1:length(replicant)
    filename=sprintf(filename_string,sample_id,pos_list(pos),repl);
    filepath=strcat(nirspectra,filename);
    import_file=importdata(filepath);
    Data_infile=import_file.data;
    wavel1=Data_infile(startrow1:endrow1,1);
    abs1=Data_infile(startrow1:endrow1,2);
    All_abs=[All_abs abs1];
 end                                                                          % end of (no-3) for loop. 
 end                                                                          % end of (no-2) for loop.
 end

%% AVERAGE CALCULATION 

Compilation_op=[];
%avgno=spectrano_sample;
%for i=1:size(Refdata,1)
  %  k=((i-1)*avgno+1);
   % av_abs=All_abs(:,k:k+(avgno-1));
   % mean_abs=mean(av_abs,2);
    Compilation_op=[Compilation_op All_abs];
%end


 %load data_matrix_for_caffeine.mat;
Wolf_Fitness_Final_ncomp = [];
for n_comp = 1:20
    
    %% Problem Definiton

    CostFunction = @(x,y,z) pls1234_modified(x,y,z);  % Cost Function

    nVar = 1;        % Dimensions of variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin =1;	% Lower Bound of Decision Variables
    VarMax = 459;    % Upper Bound of Decision Variables
    ncomp=n_comp;


%Parameters of GWO and the component of PLS
MaxIterations =20; % maximum number of iteration

PopulationSize=100;   % population of Grey Wolves
 
a =2;                %parameter a
A = 2*a*rand()-a;    %parameter A
C=2*rand();          %parameter C
%initialisation

X=zeros(MaxIterations,PopulationSize);
window_width=25;
p_value=[];
 Wolf_Fitness=[];
 Wolf_Fitness_All=[];
%CALCULATING THE FITNESS FUNCTION
for j=1:PopulationSize
    start_p=25;
    end_p=434;
    p=round((end_p-start_p)*rand()+start_p);  %FOR CHOOSING NUMBER BETWEEN 55 TO 404
    p_value=[p_value;p];     %storing the values of randomly generated row for each wolf in an array
    Start_Row_Fitness=p-24; %the start row for matrix on which pls will run
    End_Row_Fitness=p+25;   %the end row for matrix on which pls will run
    array_data=Compilation_op(Start_Row_Fitness:End_Row_Fitness,:); %the matrix on which pls will run
    Wolf_Fitness=CostFunction(array_data,Refdata,ncomp);  % the RMSECV value calculated for each such matrix or the fitness function for each wolf
    Wolf_Fitness_All=[Wolf_Fitness_All;Wolf_Fitness];  %matrix storing fittness function of all wolves whose position is generated randomly
end
    
  tblA=table(Wolf_Fitness_All,p_value);
  tblB=sortrows(tblA);
  x=table2cell(tblB);
  Wolf_Fitness_Order=cell2mat(x);
  Wolf_Fitness_Order=Wolf_Fitness_Order';
  
   %Generating position for Alpha Beta and Delta Wolves
   X_alpha=Wolf_Fitness_Order(2,1);
   X_beta=Wolf_Fitness_Order(2,2);
   X_delta=Wolf_Fitness_Order(2,3);
    
   %storing the positions generated after the random wolf position in an array
    for i=1:PopulationSize
        X(1,i)=X(1,i)+ Wolf_Fitness_Order(2,i); 
    end
    
    
    Wolf_Fitness_loop=[];
    rmsecv_component_all=[];    

%main loop of GWO
for iterations=1:MaxIterations
    array_data_loop=[];
    
    for i=1:(PopulationSize)
        D_alpha=abs(2*rand()*X_alpha-X(iterations,i));
        D_beta=abs(2*rand()*X_beta - X(iterations,i));
        D_delta=abs(2*rand()*X_delta - X(iterations,i));
        X1=X_alpha - (2*a*rand()-a)*D_alpha;
        X2=X_beta - (2*a*rand()-a)*D_beta;
        X3=X_delta - (2*a*rand()-a)*D_delta;
        X_calc=round((X1+X2+X3)/3);
        
      % DECISION FOR WOLF POSITION CHANGE OR NOT 
        
        if X_calc<25
            X(iterations+1,i)=Wolf_Fitness_Order(2,i);
            else if X_calc>434
                X(iterations+1,i)=Wolf_Fitness_Order(2,i);
                else 
                   X(iterations+1,i)=X_calc;
                end
        end          
    end
    %UPDATE FITNESS IN EACH ITERATIONS
    P_Value=[];
    Wolf_Fitness_All_loop=[];
    
    for j=1:(PopulationSize)
        p=X(iterations+1,j);
                 Start_Row_Fitness=p-24;
                 End_Row_Fitness=p+25;
                  P_Value=[P_Value;p];
         array_data_loop=Compilation_op(Start_Row_Fitness:End_Row_Fitness,:);
        
         Wolf_Fitness_loop=CostFunction(array_data_loop,Refdata,iterations);
         Wolf_Fitness_All_loop=[Wolf_Fitness_All_loop;Wolf_Fitness_loop];
         
    end
     rmsecv_component=(Wolf_Fitness_All_loop);
     rmsecv_component_all=[rmsecv_component rmsecv_component_all];
     tblC=table(Wolf_Fitness_All_loop,P_Value);     
     Wolf_Fitness_Loop_Order=[];
     tblD=sortrows(tblC);
     y=table2cell(tblD);
     Wolf_Fitness_loop_Order=cell2mat(y);
     Wolf_Fitness_loop_Order=Wolf_Fitness_loop_Order';
     
     for i=1:PopulationSize
         if Wolf_Fitness_Order(1,i)>Wolf_Fitness_loop_Order(1,i);
             Wolf_Fitness_Order(:,i)=Wolf_Fitness_loop_Order(:,i);
         
         end
     end
     tblE=table(Wolf_Fitness_Order);
     tblF=sortrows(tblE);
     z=table2cell(tblF);
     Wolf_Fitness_Final_Order=cell2mat(z);
     X_alpha=Wolf_Fitness_Final_Order(2,1);
     X_beta=Wolf_Fitness_Final_Order(2,2);
     X_delta=Wolf_Fitness_Final_Order(2,3);
    p_with_wavelength=[1:1:459;900:1.75:1701.5];                           
    clm1=X_alpha - window_width;
    predicted_start_wavel=p_with_wavelength(2,clm1);                       %%% START WAVELENGTH PREDICTION
    clm2=X_alpha + window_width;
    predicted_end_wavel=p_with_wavelength(2,clm2);
end
    Wolf_Fitness_Final_ncomp = [Wolf_Fitness_Final_ncomp;Wolf_Fitness_Final_Order];

end
