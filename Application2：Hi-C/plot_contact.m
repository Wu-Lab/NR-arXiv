function []=plot_contact(contact_mat,i,Title)

% qminHiC = quantile(contact_mat(contact_mat>0),0.01);
% qmaxHiC = quantile(contact_mat(contact_mat>0),0.99);
% 
% contact_mat(contact_mat < qminHiC) = 0;
% contact_mat(contact_mat > qmaxHiC) = qmaxHiC;

figure(i)
h=heatmap(contact_mat);
h.GridVisible = 'off';
h.Colormap =[linspace(255,255,100)' linspace(245,0,100)' linspace(245,0,100)']/255; 
h.Title=Title;
end
