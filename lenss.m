function p = vb_cw(T,S) % local function calculating brine volume of gas-free sea ice from its temperature and salinity from Cox and Weeks (1983)
    rho_i = 916.8 - 0.1403*T; % pure ice density, Pounder (1965)
    F1 = -4.732-22.45*T - 0.6397*T.^2 - 0.01074*T.^3; % Cox and Weeks (1983)
    F1(T>-2) = -4.1221*10^-2 + -1.8407*10^1*T(T>-2).^1 + 5.8402*10^-1*T(T>-2).^2 + 2.1454*10^-1*T(T>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2 = 8.903*10^-2 - 1.763*10^-2*T - 5.33*10^-4*T.^2 - 8.801*10^-6*T.^3; % Cox and Weeks (1983)
    F2(T>-2) = 9.0312*10^-2 + -1.6111*10^-2*T(T>-2).^1 + 1.2291*10^-4*T(T>-2).^2 + 1.3603*10^-4*T(T>-2).^3; % F2 from Lepparanta and Manninen (1988)
    rho_si_cox = 1./((F1 - rho_i.*S/1000.*F2)./(rho_i.*F1)); % Cox and Weeks (1983)
    p = rho_si_cox .* S ./ F1 / 1000; % Cox and Weeks (1983)
end

close all; clc; clear; load('lenss_st6.mat');
t(1) = datetime('14-Jan-2020'); lon = 11.4604; lat = -69.6629; % Station 6 time and location
Srho = gsw_SA_from_SP(SPrho,0,lon,lat); S = gsw_SA_from_SP(SP,0,lon,lat); % absolute/practical salinity coversion from TEOS-10
clearvars -except hS hT hrho dS dT drho zT zzS zzrho T T_lab S Srho rho t event
zS = mean(zzS,2) * hS / zzS(end,2) / 100;
T_S = interp1(zT * hT / zT(end),T,zS*100,'linear','extrap');
rho_i = 916.8 - 0.1403*T_S; % pure ice density, Pounder (1965)
F1 = -4.732-22.45*T_S - 0.6397*T_S.^2 - 0.01074*T_S.^3; % Cox and Weeks (1983)
F2 = 8.903*10^-2 - 1.763*10^-2*T_S - 5.33*10^-4*T_S.^2 - 8.801*10^-6*T_S.^3; % Cox and Weeks (1983)
F1(T_S>-2) = -4.1221*10^-2 + -1.8407*10^1*T_S(T_S>-2).^1 + 5.8402*10^-1*T_S(T_S>-2).^2 + 2.1454*10^-1*T_S(T_S>-2).^3; % F1 from Lepparanta and Manninen (1988)
F2(T_S>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_S(T_S>-2).^1 + 1.2291*10^-4*T_S(T_S>-2).^2 + 1.3603*10^-4*T_S(T_S>-2).^3; % F2 from Lepparanta and Manninen (1988)
rho_si_cox = 1./((F1 - rho_i.*S/1000.*F2)./(rho_i.*F1)); % Cox and Weeks (1983)
vb_S_cw = rho_si_cox .* S ./ F1; vb_S_cw = vb_S_cw/1000; vb_S_cw(vb_S_cw > 0.5) = NaN; vb_S_cw(vb_S_cw < 0) = NaN; vb_S = vb_S_cw; % Cox and Weeks (1983)
Sb = -1.2 - 21.8*T_S - 0.919*T_S.^2 - 0.0178*T_S.^3; % brine salinity for seawater (Assur, 1958)
Sb_van = -18.7*T_S - 0.519*T_S.^2 - 0.00535*T_S.^3; % brine salinity for seawater (Vancopenolle, 2019, eq. 10)
rho_w = 1000.3 + 0.78237 * Sb + 2.8008*10^-4 * Sb.^2; % seawater density (Schwerdtfeger, 1963) 
rho_b = 1000.3 + 0.78237 * Sb_van + 2.8008*10^-4 * Sb_van.^2; % seawater density (Schwerdtfeger, 1963) 
vb_S_van = ( 1 + (Sb_van./S - 1) .* rho_b./rho_i).^-1; % brine volume for gas-free ice without minerals (Vancopenolle, 2019, eq. 6)
phi_v = (1 - S./Sb)./(1+S./Sb.*(rho_i./rho_w-1)); vb_S_notz = 1 - phi_v; vb_S_notz(vb_S_notz > 0.6) = NaN; vb_S_notz(vb_S_notz < 0) = NaN; % volume solid fraction (Notz, 2005)
zrho = mean(zzrho,2) * hrho / zzrho(end,2) / 100;
T_rho = interp1(zT * hT / zT(end),T,zrho*100,'linear','extrap');
F1_pr_rho = -4.732-22.45*T_lab - 0.6397*T_lab.^2 - 0.01074*T_lab.^3;
F2_pr_rho = 8.903*10^-2 - 1.763*10^-2*T_lab - 5.33*10^-4*T_lab.^2 - 8.801*10^-6*T_lab.^3;
vb_pr_rho = rho .* Srho ./ F1_pr_rho; % brine volume for T_lab
rhoi_pr = (917-1.403*10^-1*T_lab); % pure ice density, Pounder (1965)
vg_pr = max((1-rho.*(F1_pr_rho-rhoi_pr.*Srho/1000.*F2_pr_rho)./(rhoi_pr.*F1_pr_rho)),0,'includenan'); % gas volume for T_lab
F3_pr = rhoi_pr.*Srho/1000./(F1_pr_rho-rhoi_pr.*Srho/1000.*F2_pr_rho);
rhoi_rho = (917-1.403*10^-1*T_rho); % pure ice density, Pounder (1965) for T_insitu
F1_rho = -4.732-22.45*T_rho - 0.6397*T_rho.^2 - 0.01074*T_rho.^3; % Cox and Weeks (1983)
F1_rho(T_rho>-2) = -4.1221*10^-2 + -1.8407*10^1*T_rho(T_rho>-2).^1 + 5.8402*10^-1*T_rho(T_rho>-2).^2 + 2.1454*10^-1*T_rho(T_rho>-2).^3; % F1 from Lepparanta and Manninen (1988)
F2_rho = 8.903*10^-2 - 1.763*10^-2*T_rho - 5.33*10^-4*T_rho.^2 - 8.801*10^-6*T_rho.^3; % Cox and Weeks (1983)
F2_rho(T_rho>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_rho(T_rho>-2).^1 + 1.2291*10^-4*T_rho(T_rho>-2).^2 + 1.3603*10^-4*T_rho(T_rho>-2).^3; % F2 from Lepparanta and Manninen (1988)
F3_rho = rhoi_rho.*Srho/1000./(F1_rho-rhoi_rho.*Srho/1000.*F2_rho);
vb_rho = vb_pr_rho .* F1_pr_rho ./ F1_rho / 1000; vb_rho(vb_rho > 0.6) = NaN; vb_rho(vb_rho < 0) = NaN;  % Brine volume for T_insitu, CW + LM
vg = max(0,(1-(1-vg_pr).*(rhoi_rho./rhoi_pr).*(F3_pr.*F1_pr_rho./F3_rho./F1_rho)),'includenan'); % Gas volume for T_insitu, CW + LM
rho_si = (-vg_pr+1).*rhoi_rho.*F1_rho./(F1_rho-rhoi_rho.*Srho/1000.*F2_rho); rho_si(isnan(vb_rho)) = NaN; rho_si(isnan(vb_rho)) = NaN; % density
rho_si_gf = (1).*rhoi_rho.*F1_rho./(F1_rho-rhoi_rho.*Srho/1000.*F2_rho); rho_si(isnan(vb_rho)) = NaN; rho_si(isnan(vb_rho)) = NaN; % density
% clearvars -except hS hT hrho dS dT drho zT zS zzS zzrho zrho T T_lab T_S S Srho rho rho_si vb_S vb_rho vg vg_pr t event vb_S_van rho_si_gf

figure
tile = tiledlayout(1,5); tile.TileSpacing = 'compact'; tile.Padding = 'none';
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; % colors

nexttile
plot(S,zS,Srho,zrho);
set(gca, 'YDir','reverse'); leg = legend('SAL','DEN','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.3,18*0.3];
hXLabel = xlabel('Salinity (g kg^-^1)'); hYLabel = ylabel('Ice thickness (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
p = text(8,1.4,sprintf('%.1f ± %.1f',mean(S),std(S))); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8);
p = text(8,1.5,sprintf('%.1f ± %.1f',mean(Srho),std(Srho))); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',8);

nexttile
plot(T,zT/100);
set(gca, 'YDir','reverse');
hXLabel = xlabel('Temperature (°C)'); set([hXLabel gca],'FontSize',8,'FontWeight','normal'); ylim([0 2]);
p = text(-2.45,1.45,sprintf('%.1f ± %.1f',mean(T),std(T))); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',8);

nexttile
plot(rho_si,zrho,rho,zrho);
set(gca, 'YDir','reverse'); leg = legend('T','Tlab','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.3,18*0.3];
hXLabel = xlabel('Density (kg m^-^3)'); set([hXLabel gca],'FontSize',8,'FontWeight','normal');
p = text(855,1.4,sprintf('%.0f ± %.0f',mean(rho_si),std(rho_si))); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',8);
p = text(855,1.5,sprintf('%.0f ± %.0f',mean(rho),std(rho))); set(p,'Color',c{2},'HorizontalAlignment','left','FontSize',8);

nexttile
plot(vb_S,zS,vb_rho,zrho); hold on
% plot(vb_cw(T_S,S)+0.01,zS,'k--'); % using vb_cw = vb_cw(T,S);
% set(gca, 'YDir','reverse'); leg = legend('SAL, CW','SAL, VC','DEN','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.3,18*0.3];
set(gca, 'YDir','reverse'); leg = legend('SAL','DEN','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.3,18*0.3];
hXLabel = xlabel('Brine volume'); set([hXLabel gca],'FontSize',8,'FontWeight','normal');
p = text(0.19,1.4,sprintf('%.2f ± %.2f',mean(vb_S),std(vb_S))); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8);
p = text(0.19,1.5,sprintf('%.2f ± %.2f',mean(vb_rho),std(vb_rho))); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',8);

nexttile
plot(vg,zrho,vg_pr,zrho);
set(gca, 'YDir','reverse'); leg = legend('T','Tlab','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.3,18*0.3];
hXLabel = xlabel('Gas volume'); set([hXLabel gca],'FontSize',8,'FontWeight','normal');
p = text(0.095,1.4,sprintf('%.2f ± %.2f',mean(vg),std(vg))); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8);
p = text(0.095,1.5,sprintf('%.2f ± %.2f',mean(vg_pr),std(vg_pr))); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',8);

annotation('textbox',[0.005 .51 0.01 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.144 .51 0.15 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.273 .51 0.28 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.406 .51 0.41 .51],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.534 .51 0.54 .51],'String','(e)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');

clearvars tile leg hXLabel hYLabel