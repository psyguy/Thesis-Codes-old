load_system('BiNode');
hws = get_param('BiNode', 'modelworkspace');

P1 = 7;
P2 = 4.5;
%{
PFC: P1=4.5, P2 = 0;
PAC: P1=7, P2 = 0;
AAC: P1=7, P2 = 7;
FFC: P1=4.5, P2 = 4.5;
AFC: P1=4.5, P2 = 7;

%}
tic;
Set_ModelParam(hws,'P1',P1);
Set_ModelParam(hws,'P2',P2);
sim('BiNode');
ResultAnalysis;
toc;