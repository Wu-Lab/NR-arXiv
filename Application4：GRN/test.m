%% NR-B improves the accuracy of Gene Regulatory Networks inference 

clear
clc
addpath(genpath(pwd));

%%
auc_before=zeros(1,9);
auc_after_ND=zeros(1,9);
auc_after_NR=zeros(1,9);
pr_before=zeros(1,9);
pr_after_ND=zeros(1,9);
pr_after_NR=zeros(1,9);
methods={'CLR';'ARACNE';'MI';'Pearson';'Spearman';'GENIE3';'TIGRESS';'Inferelator';'ANOVA'};
m=2;

disp('******************* Evalution: *******************')

for j=1:9
    
    disp(['method name:',cell2mat(methods(j))])
    
    %% load data
    input_network=load_data(1,cell2mat(methods(j)));
    gold_file = 'networks/DREAM5_NetworkInference_Edges_Network1.tsv';
    gold_edges = load_dream_network(gold_file);
    
    %% compute performance before NR
    prediction_before=change_network_format(input_network);
    [auc_before(j), pr_before(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_before);
    
    %% apply NR algorithm and ND algorithm
    input_network = (input_network-min(min(input_network)))./(max(max(input_network))-min(min(input_network)));
    output_network_ND = ND_regulatory(input_network);
    output_network_NR = NR_B(input_network,m);
    
    %% compute performance after NR
    prediction_after_ND=change_network_format(output_network_ND);
    [auc_after_ND(j), pr_after_ND(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_ND);
    
    prediction_after_NR=change_network_format(output_network_NR);
    [auc_after_NR(j), pr_after_NR(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_NR);
    
end



auc_score=[auc_before;auc_after_ND;auc_after_NR];
pr_score=[pr_before;pr_after_ND;pr_after_NR];


set(gcf,'color','white')
figure(1)
h=bar(auc_score');
legend('original','ND','NR-B')
set(h(1),'facecolor','#999999')
set(h(2),'facecolor','#009E73')
set(h(3),'facecolor','#F0E442')
ylabel('AUROC')
ylim([0,1])
set(gca,'XTickLabel',methods)
set(gca,'XTickLabelRotation',90)
set(gca,'FontSize',12,'FontName','bold')
set(h,'edgecolor','none');
grid on;


set(gcf,'color','white')
figure(2)
h=bar(pr_score');
legend('original','ND','NR-B')
set(h(1),'facecolor','#999999')
set(h(2),'facecolor','#009E73')
set(h(3),'facecolor','#F0E442')
ylabel('AUPR')
set(gca,'XTickLabel',methods)
set(gca,'XTickLabelRotation',90)
set(gca,'FontSize',12,'FontName','bold')
set(h,'edgecolor','none');
grid on;

