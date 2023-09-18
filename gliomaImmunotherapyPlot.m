%-------------------------------------------------------------------%
%--------------------------- Plot Solution -------------------------%
%-------------------------------------------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;
for i=1:length(output.meshhistory);
  mesh(i).meshPoints = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).time =  output.meshhistory(i).result.solution.phase.time;
  mesh(i).iteration = i*ones(size(mesh(i).meshPoints));
  mesh(i).iterationTime = i*ones(size(mesh(i).time));
end;

tf = solution.phase.time(end);

%Legendlabel = ["Cancer"; "T cells"; "MDSCs"];
fs = 18; % font size
ms = 8; % marker size
lw = 1.25; % line width

finaltime = t0Max + bounds.phase.duration.upper;
yaxisub = [5e5, 5e7, 5e4, 1, 1]; % UB for y-axis


figure(1)
stateeqns = tiledlayout(2,3,'TileSpacing','compact');
%ax = nexttile(1);
nexttile(1)
% cancer = plot(solution.phase(1).time,solution.phase(1).state(:,1),'-o');
% hold on
% tcell = plot(solution.phase(1).time,solution.phase(1).state(:,2),'-o');
% mdsc = plot(solution.phase(1).time,solution.phase(1).state(:,3),'-o');
% hold off
pp = plot(solution.phase(1).time,solution.phase(1).state(:,1),'-o');
%axis([0 finaltime 0 yaxisub(1)])
xl = xlabel('$t$','Interpreter','LaTeX');
%xl = xlabel('$q(t)$','Interpreter','LaTeX');
yl = ylabel('C(t)','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
%set(cancer,'LineWidth',lw,'MarkerSize',ms);
%set(tcell,'LineWidth',lw,'MarkerSize',ms);
%set(mdsc,'LineWidth',lw,'MarkerSize',ms);
%grid on;
%print("gliomaImmunotherapyState1.eps","-depsc")
%print -depsc2 gliomaImmunotherapyState1.eps
%exportgraphics(gcf,"gliomaImmunotherapyState1.eps","ContentType","vector")

nexttile(2)
pp = plot(solution.phase(1).time,solution.phase(1).state(:,2),'-o');
%axis([0 finaltime 0 yaxisub(2)])
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
% %print -depsc2 gliomaImmunotherapyState2.eps
% %saveas(gcf, 'gliomaImmunotherapyState2.eps')

nexttile(3)
pp = plot(solution.phase(1).time,solution.phase(1).state(:,3),'-o');
%axis([0 finaltime 0 yaxisub(3)])
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$M(t)$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
% %print -depsc2 gliomaImmunotherapyState3.eps

nexttile(4)
pp = plot(solution.phase(1).time,solution.phase(1).control(:,1),'-o');
%axis([0 finaltime 0 yaxisub(4)])
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\gamma u_1$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
%print -depsc2 gliomaImmunotherapyControl1.eps

nexttile(5)
pp = plot(solution.phase(1).time,solution.phase(1).control(:,2),'-o');
%axis([0 finaltime 0 yaxisub(5)])
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
% %print -depsc2 gliomaImmunotherapyControl2.eps

set(gcf,'Position',[0 300 1000 650],'PaperPositionMode','auto');
%lg = legend(ax,Legendlabel,'FontSize',fs);
%lg.Layout.Tile = 'east';

%% Costate Equations (aka Adjoint Equations)

figure(2)
adjointeqns = tiledlayout(1,3,'TileSpacing','compact');
nexttile(1)
pp = plot(solution.phase(1).time,solution.phase(1).costate(:,1),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_C(t)$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
%print -depsc2 gliomaImmunotherapyCostate1.eps


nexttile(2)
pp = plot(solution.phase(1).time,solution.phase(1).costate(:,2),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
%print -depsc2 gliomaImmunotherapyCostate2.eps

nexttile(3)
pp = plot(solution.phase(1).time,solution.phase(1).costate(:,3),'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_M(t)$','Interpreter','LaTeX');
set(xl,'FontSize',fs);
set(yl,'FontSize',fs);
set(gca,'FontSize',fs-2,'FontName','Times');
set(pp,'LineWidth',lw,'MarkerSize',ms);
grid on;
%print -depsc2 gliomaImmunotherapyCostate3.eps
set(gcf,'Position',[0 300 900 350],'PaperPositionMode','auto');

%% Mesh

figure(3)
%mesh = tiledlayout(1,2,'TileSpacing','compact');
%nexttile(1)
for i=1:length(mesh);
  pp = plot(mesh(i).meshPoints*tf,mesh(i).iteration,'bo');
  set(pp,'LineWidth',lw);
  hold on;
end;
xl = xlabel('Mesh Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',fs);
set(yl,'Fontsize',fs);
set(gca,'YTick',0:length(mesh),'FontSize',fs-2);
grid on;
%print -depsc2 gliomaImmunotherapyMeshRefinement.eps

figure(4)
%nexttile(2)
for i=1:length(mesh);
  pp = plot(mesh(i).time,mesh(i).iterationTime,'bo');
  set(pp,'LineWidth',lw);
  hold on;
end;
xl = xlabel('Collocation Point Location');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',fs);
set(yl,'Fontsize',fs);
set(gca,'YTick',0:length(mesh),'FontSize',fs-2);
grid on;
%print -depsc2 gliomaImmunotherapyMeshRefinementTime.eps