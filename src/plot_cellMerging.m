close all
clc
clear all

% Get screen size
scrsz = get(0,'screensize');
scrsz2 = scrsz;
scrsz2(4) = scrsz(4)/2; % square
scrsz2(3) = scrsz2(4);

nelem = [8,16,32,64,128];
dx = 1./nelem; 
ndof = zeros(5,5); 
for i = 1:5
    ndof(i,:) = nelem*(i+1);
end

%% m2start = -0.5 - dx(2)*0.99
% Baseline Abutting Method
base.legendre.p1.l2 = [1.6192323442766485E-002   5.5463743682936896E-003   1.1523940541828092E-003   1.9492368358601888E-004   3.4677255956864258E-005];
base.legendre.p1.cons = [2.1402403156450811E-004  -5.5089999743261620E-005  -5.5814355262792659E-006  -2.5117508424543633E-007  -3.8742272931435728E-009];

base.legendre.p2.l2 = [3.1172035890003353E-003   2.5836722634559685E-004   2.0358670146593475E-005   2.3144704192818786E-006   2.8616992403088300E-007];
base.legendre.p2.cons = [-1.2732668244837264E-004  -1.5343553475613558E-006  -2.1277249420692534E-008  -7.7716750396250589E-010  -4.4681973399018915E-011];

base.legendre.p3.l2 = [4.6738078468553133E-004   1.6852588106675376E-005   8.4821948988185273E-007   5.2238918471511937E-008   1.2756984616776392E-008];
base.legendre.p3.cons = [1.0270445788065907E-006   7.3374229053735363E-009  -6.4898900142829063E-010   7.4165702745876416E-010   8.0106703909121890E-010];

base.legendre.p4.l2 = [4.9171322299930194E-005   1.1180711012635708E-006   3.8098986030597548E-008   1.3357828582365599E-008   1.1566783257513384E-008];
base.legendre.p4.cons = [5.1535640198385746E-007   2.4886523605771593E-009   3.4894656192330231E-009   3.5039834089811706E-009   3.5042774515492425E-009];

base.legendre.p5.l2 = [6.7157097371632029E-006   9.5844108831121343E-008   1.9641317833152333E-008   1.3931055256015576E-008   1.2526408489692476E-008];
base.legendre.p5.cons = [-3.4980001256679216E-009   2.9228249334711798E-009   2.8813135627414610E-009   2.8810332244888492E-009   2.8802265225613688E-009];

base.lagrange.p1.l2 = [0.13905677101460567        5.2522707400047830E-002   1.1371263253338055E-002   1.9465140397926721E-003   3.4709809909961377E-004];
base.lagrange.p1.cons = [2.5809996556810377E-003  -5.4506782951047050E-004  -6.2567924626277804E-005  -3.6065007995578213E-006  -1.8312580379875598E-007];

base.lagrange.p2.l2 = [2.8715542958208582E-002   2.6916331092555821E-003   2.1234737319991985E-004   2.4069049414627405E-005   2.9698656840838573E-006];
base.lagrange.p2.cons = [-1.3826073488102386E-003  -2.8014648676355502E-005  -1.2707673529566144E-006  -8.2669196599383099E-008  -5.3819386169884353E-009];

base.lagrange.p3.l2 = [4.7190705848199109E-003   1.7618488181622786E-004   8.9339422387308516E-006   5.6378424937694400E-007   1.2319245960449388E-007];
base.lagrange.p3.cons = [-9.7928305631267420E-07  -5.8863419977184250E-007  -3.0490229674384750E-008   6.6429308720472591E-009   7.9865597468753435E-009];

base.lagrange.p4.l2 = [5.1546559875049211E-004   1.1467957216868954E-005   3.7857455005431421E-007   1.8352725989409218E-007   1.1121466206592967E-007];
base.lagrange.p4.cons = [3.7249135267369837E-006  -1.5509117878842460E-008   3.4038895424526316E-008   3.5025322531190284E-008   3.5042550861064115E-008];

base.lagrange.p5.l2 = [7.2583458104059359E-005   1.0357192468086048E-006   1.4008420017077255E-007   1.7434504909779441E-007   1.4993191502827124E-007];
base.lagrange.p5.cons = [-2.3261550174957790E-007   2.8015271491543814E-008   2.8804939655735495E-008   2.8810566876646959E-008   2.8802266516247954E-008];

% Conservative Overset WITH cell merging

cons.legendre.p1.l2 = [1.5343867248100811E-002   4.9527783109428562E-003   9.1511641008077576E-004   1.2837168104084824E-004   1.7950048525372715E-005];
cons.legendre.p1.cons = [-9.1401686055925246E-007  -8.8376291847158672E-008  -4.2944623968033113E-009   9.5057671456455495E-010   8.1442088489813713E-010];

cons.legendre.p2.l2 = [3.0856774565902233E-003   2.0983997033755593E-004   8.6303908350958463E-006   5.8468668055943048E-007   1.5344243887565960E-007];
cons.legendre.p2.cons = [-5.4618390643579318E-008  -5.7510571235819263E-009  -3.1144142126349195E-009  -2.1074259767828174E-009  -1.2337492277803008E-009];

cons.legendre.p3.l2 = [3.5317521845952771E-004   7.5465162545232008E-006   3.7182215975436738E-007   2.1007049446592616E-007   1.7390125794548207E-007];
cons.legendre.p3.cons = [-1.2349907894948497E-007  -8.4690109092289134E-011   3.6569894265592318E-009   2.6092678989675377E-009   1.5316150986377508E-009];

cons.legendre.p4.l2 = [3.4263040992898624E-005   3.7061320813902521E-007   2.4491141395242653E-007   2.3908413823671063E-007   1.9822081502381796E-007];
cons.legendre.p4.cons = [-3.3117729392184714E-010  -4.6735204456416213E-009  -4.3054433421230165E-009  -2.9763325623544112E-009  -1.7447722364116736E-009];

cons.legendre.p5.l2 = [3.3725834680316162E-006   2.0117279770511893E-007   2.6478482271882418E-007   2.5862698957822409E-007   2.1440616753902815E-007];
cons.legendre.p5.cons = [3.6081780237229300E-008   4.8850695918978282E-009   4.6397280392684870E-009   3.2086964751876224E-009   1.8805727669546357E-009];

cons.lagrange.p1.l2 = [0.13091100162057964        4.6743460591335113E-002   9.0282272507577113E-003   1.2841419885300939E-003   1.8039929898380172E-004];
cons.lagrange.p1.cons = [-2.9868863937743129E-007  -3.1733510319487124E-008  -3.7768582705766107E-009  -4.6558240596183964E-010  -5.7957860732926747E-011];

cons.lagrange.p2.l2 = [2.8166788897167321E-002   2.2296103048719271E-003   1.0501062405271639E-004   8.6378569183109159E-006   1.0179298476077247E-006];
cons.lagrange.p2.cons = [4.4095609719185802E-008   1.4714956475536667E-009   4.9208914720821895E-011   1.7256751583261121E-012   6.5336624999190462E-014];

cons.lagrange.p3.l2 = [3.5737165469216766E-003   9.1198535593244387E-005   4.1745302116984371E-006   3.0387014059787933E-007   1.1964167225031556E-007];
cons.lagrange.p3.cons = [-1.1853670822681295E-008  -4.0376163523703212E-010  -1.2790268844042885E-011  -4.0056846728475648E-013  -1.1823875212257917E-014];

cons.lagrange.p4.l2 = [3.7578277688742913E-004   4.1010088672093868E-006   1.7903703576211899E-007   1.7767452689355275E-007   1.0106540154787511E-007];
cons.lagrange.p4.cons = [-6.0015614700148490E-010  -3.1335489758532731E-012  -1.4988010832439613E-014   2.2204460492503131E-016   7.2164496600635175E-016];

cons.lagrange.p5.l2 = [4.2483859342141195E-005   5.8105031583470969E-007   1.3478763641681851E-007   1.7013187189705483E-007   1.4489828341683108E-007];
cons.lagrange.p5.cons = [1.5315521073588911E-010   7.6932904491400222E-013   5.4400928206632670E-015  -2.7755575615628914E-016  -5.5511151231257827E-017];

%% Plots

% % figure(1)
% % subplot(1,2,1);
% % loglog(ndof(1,:),abs(base.lagrange.p1.l2),'b--',ndof(1,:),cons.lagrange.p1.l2,'b-','linewidth',3)
% % hold on; 
% % loglog(ndof(2,:),abs(base.lagrange.p2.l2),'r--',ndof(2,:),cons.lagrange.p2.l2,'r-','linewidth',3)
% % loglog(ndof(3,:),abs(base.lagrange.p3.l2),'g--',ndof(3,:),cons.lagrange.p3.l2,'g-','linewidth',3)
% % loglog(ndof(4,:),abs(base.lagrange.p4.l2),'m--',ndof(4,:),cons.lagrange.p4.l2,'m-','linewidth',3); 
% % loglog(ndof(5,:),abs(base.lagrange.p5.l2),'c--',ndof(5,:),cons.lagrange.p5.l2,'c-','linewidth',3); 
% % % loglog(dx,(base.lagrange.p1.l2(2)/dx(2).^2)*dx.^2,'k--')
% % % loglog(dx,(base.lagrange.p2.l2(2)/dx(2).^3)*dx.^3,'k--')
% % % loglog(dx,(base.lagrange.p3.l2(2)/dx(2).^4)*dx.^4,'k--')
% % % loglog(dx,(base.lagrange.p4.l2(2)/dx(2).^5)*dx.^5,'k--')
% % % loglog(dx,(base.lagrange.p5.l2(2)/dx(2).^6)*dx.^6,'k--')
% % grid on;
% % legend('base p=1','cons p=1','base p=2','cons p=2','base p=3','cons p=3','base p=4','cons p=4', 'base p=5','cons p=5','location','southeast');
% % grid minor; 
% % xlabel('1/N_{elem}','fontsize',18)
% % ylabel('L2 Error','fontsize',18)
% % set(gca,'fontsize',18')
% % % ylim([1e-8 1])
% % title({'Lagrange','99% Cut'},'FontSize',30')
% % 
% % subplot(1,2,2)
% % loglog(ndof(1,:),abs(base.lagrange.p1.cons),'b--',ndof(1,:),abs(cons.lagrange.p1.cons),'b-','linewidth',3)
% % hold on; 
% % loglog(ndof(2,:),abs(base.lagrange.p2.cons),'r--',ndof(2,:),abs(cons.lagrange.p2.cons),'r-','linewidth',3)
% % loglog(ndof(3,:),abs(base.lagrange.p3.cons),'g--',ndof(3,:),abs(cons.lagrange.p3.cons),'g-','linewidth',3)
% % loglog(ndof(4,:),abs(base.lagrange.p4.cons),'m--',ndof(4,:),abs(cons.lagrange.p4.cons),'m-','linewidth',3); 
% % loglog(ndof(5,:),abs(base.lagrange.p5.cons),'c--',ndof(5,:),abs(cons.lagrange.p5.cons),'c-','linewidth',3); 
% % % loglog(dx,abs(base.lagrange.p1.cons(2)/dx(2).^2)*dx.^2,'k--')
% % % loglog(dx,abs(base.lagrange.p2.cons(3)/dx(3).^3)*dx.^3,'k--')
% % % loglog(dx,abs(base.lagrange.p3.cons(3)/dx(3).^4)*dx.^4,'k--')
% % % loglog(dx,abs(base.lagrange.p4.cons(2)/dx(2).^5)*dx.^5,'k--')
% % % loglog(dx,abs(base.lagrange.p5.cons(2)/dx(2).^6)*dx.^6,'k--')
% % grid on;
% % grid minor; 
% % xlabel('1/N_{elem}','fontsize',18)
% % ylabel('Conservation Error','fontsize',18)
% % set(gca,'fontsize',18')
% % title({'Lagrange','99% Cut'},'FontSize',30')


figure(10)
subplot(1,2,1);
loglog(dx,abs(base.lagrange.p1.l2),'b--',dx,cons.lagrange.p1.l2,'b-','linewidth',3)
hold on; 
loglog(dx,abs(base.lagrange.p2.l2),'r--',dx,cons.lagrange.p2.l2,'r-','linewidth',3)
loglog(dx,abs(base.lagrange.p3.l2),'g--',dx,cons.lagrange.p3.l2,'g-','linewidth',3)
loglog(dx,abs(base.lagrange.p4.l2),'m--',dx,cons.lagrange.p4.l2,'m-','linewidth',3); 
loglog(dx,abs(base.lagrange.p5.l2),'c--',dx,cons.lagrange.p5.l2,'c-','linewidth',3); 
loglog(dx,(base.lagrange.p1.l2(2)/dx(2).^2)*dx.^2,'k--')
loglog(dx,(base.lagrange.p2.l2(2)/dx(2).^3)*dx.^3,'k--')
loglog(dx,(base.lagrange.p3.l2(2)/dx(2).^4)*dx.^4,'k--')
loglog(dx,(base.lagrange.p4.l2(2)/dx(2).^5)*dx.^5,'k--')
loglog(dx,(base.lagrange.p5.l2(2)/dx(2).^6)*dx.^6,'k--')
grid on;
legend('base p=1','cons p=1','base p=2','cons p=2','base p=3','cons p=3','base p=4','cons p=4', 'base p=5','cons p=5','location','southeast');
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('L2 Error','fontsize',18)
set(gca,'fontsize',18')
ylim([1e-8 1])
title({'Lagrange','99% Cut'},'FontSize',30')

subplot(1,2,2)
loglog(dx,abs(base.lagrange.p1.cons),'b--',dx,abs(cons.lagrange.p1.cons),'b-','linewidth',3)
hold on; 
loglog(dx,abs(base.lagrange.p2.cons),'r--',dx,abs(cons.lagrange.p2.cons),'r-','linewidth',3)
loglog(dx,abs(base.lagrange.p3.cons),'g--',dx,abs(cons.lagrange.p3.cons),'g-','linewidth',3)
loglog(dx,abs(base.lagrange.p4.cons),'m--',dx,abs(cons.lagrange.p4.cons),'m-','linewidth',3); 
loglog(dx,abs(base.lagrange.p5.cons),'c--',dx,abs(cons.lagrange.p5.cons),'c-','linewidth',3); 
% loglog(dx,abs(base.lagrange.p1.cons(2)/dx(2).^2)*dx.^2,'k--')
% loglog(dx,abs(base.lagrange.p2.cons(3)/dx(3).^3)*dx.^3,'k--')
% loglog(dx,abs(base.lagrange.p3.cons(3)/dx(3).^4)*dx.^4,'k--')
% loglog(dx,abs(base.lagrange.p4.cons(2)/dx(2).^5)*dx.^5,'k--')
% loglog(dx,abs(base.lagrange.p5.cons(2)/dx(2).^6)*dx.^6,'k--')
grid on;
ylim([1e-17 1e-2])
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('Conservation Error','fontsize',18)
set(gca,'fontsize',18')
title({'Lagrange','99% Cut'},'FontSize',30')
set(gcf,'Position',scrsz);
saveas(gcf,'ConsOverset_CellMerging_Lagrange_99pctCut.png');
close(gcf);


% 
% figure(2)
% subplot(1,2,1);
% loglog(dx,abs(base.legendre.p1.l2),'b--',dx,cons.legendre.p1.l2,'b','linewidth',3)
% hold on; 
% loglog(dx,abs(base.legendre.p2.l2),'r--',dx,cons.legendre.p2.l2,'r','linewidth',3)
% loglog(dx,abs(base.legendre.p3.l2),'g--',dx,cons.legendre.p3.l2,'g','linewidth',3)
% loglog(dx,abs(base.legendre.p4.l2),'m--',dx,cons.legendre.p4.l2,'m','linewidth',3); 
% loglog(dx,abs(base.legendre.p5.l2),'c--',dx,cons.legendre.p5.l2,'c','linewidth',3); 
% loglog(dx,(base.legendre.p1.l2(2)/dx(2).^2)*dx.^2,'k--')
% loglog(dx,(base.legendre.p2.l2(2)/dx(2).^3)*dx.^3,'k--')
% loglog(dx,(base.legendre.p3.l2(2)/dx(2).^4)*dx.^4,'k--')
% loglog(dx,(base.legendre.p4.l2(2)/dx(2).^5)*dx.^5,'k--')
% loglog(dx,(base.legendre.p5.l2(2)/dx(2).^6)*dx.^6,'k--')
% grid on;
% legend('base p=1','cons p=1','base p=2','cons p=2','base p=3','cons p=3','base p=4','cons p=4', 'base p=5','cons p=5','location','southeast');
% grid minor; 
% xlabel('1/N_{elem}','fontsize',18)
% ylabel('L2 Error','fontsize',18)
% set(gca,'fontsize',18')
% % ylim([1e-8 1])
% title({'Legendre','99% Cut'},'FontSize',30')
% 
% subplot(1,2,2)
% loglog(dx,abs(base.legendre.p1.cons),'b--',dx,abs(cons.legendre.p1.cons),'b','linewidth',3)
% hold on; 
% loglog(dx,abs(base.legendre.p2.cons),'r--',dx,abs(cons.legendre.p2.cons),'r','linewidth',3)
% loglog(dx,abs(base.legendre.p3.cons),'g--',dx,abs(cons.legendre.p3.cons),'g','linewidth',3)
% loglog(dx,abs(base.legendre.p4.cons),'m--',dx,abs(cons.legendre.p4.cons),'m','linewidth',3); 
% loglog(dx,abs(base.legendre.p5.cons),'c--',dx,abs(cons.legendre.p5.cons),'c','linewidth',3); 
% loglog(dx,abs(base.legendre.p1.cons(2)/dx(2).^2)*dx.^2,'k--')
% loglog(dx,abs(base.legendre.p2.cons(3)/dx(3).^3)*dx.^3,'k--')
% loglog(dx,abs(base.legendre.p3.cons(3)/dx(3).^4)*dx.^4,'k--')
% loglog(dx,abs(base.legendre.p4.cons(2)/dx(2).^5)*dx.^5,'k--')
% loglog(dx,abs(base.legendre.p5.cons(2)/dx(2).^6)*dx.^6,'k--')
% grid on;
% grid minor; 
% xlabel('1/N_{elem}','fontsize',18)
% ylabel('Conservation Error','fontsize',18)
% set(gca,'fontsize',18')
% title({'Legendre','99% Cut'},'FontSize',30')

%% m2start = -0.5 - dx(2)*0.90
% Baseline Abutting Method
base.legendre.p1.l2 = [1.6765449067581831E-002   5.7242656668346961E-003   1.2206951905958544E-003   2.1365249958722176E-004   3.9464805396971161E-005];
base.legendre.p1.cons = [-9.4496309399845491E-005  -6.2303509111988342E-005  -2.0919789927367760E-006   9.6643524169492956E-008   2.7298075633286345E-008];

base.legendre.p2.l2 = [3.2090886436799268E-003   2.4424042267942304E-004   1.5492158106219127E-005   1.5806009481025005E-006   1.9229961664220241E-007];
base.legendre.p2.cons = [-1.8315247259234613E-004  -8.6545720316688612E-007   2.2816193587238409E-009  -4.3953660155970908E-012   6.0920296585109668E-012];

base.legendre.p3.l2 = [5.0337548393662463E-004   1.8860531301872488E-005   9.8357178104406351E-007   6.0868842247452506E-008   1.2936797534971091E-008];
base.legendre.p3.cons = [-9.5798795170798368E-006   5.6432799037942072E-008   1.8292017275878258E-009   3.4347663602218859E-011  -9.7384669772715426E-012];

base.legendre.p4.l2 = [4.4011746398404877E-005   7.0077841142936975E-007   2.5723532283141137E-008   1.2707483254631212E-008   1.0830623355647697E-008];
base.legendre.p4.cons = [-3.0264019609815840E-007  -2.7169601221621598E-009   1.6791833340468187E-009   1.7100567697991842E-009   1.7103266164442132E-009];

base.legendre.p5.l2 = [7.9469426215944998E-006   1.1532338271007949E-007   2.0025149463224594E-008   1.4437568094157017E-008   1.3079155483203415E-008];
base.legendre.p5.cons = [-1.7259153542448225E-007   4.5347092678116319E-009   4.0502249887275354E-009   4.0284327262885000E-009   4.0225528810666766E-009];

base.lagrange.p1.l2 = [0.14445001417793107        5.4256555024397628E-002   1.2049762304374897E-002   2.1334221313525480E-003   3.9489362306799016E-004];
base.lagrange.p1.cons = [-3.5090952909039297E-004  -5.9195613677087700E-004  -2.3242288216429596E-005   5.5120368464001146E-007   2.1747653428105806E-007];

base.lagrange.p2.l2 = [2.9744245998739302E-002   2.5607220579665761E-003   1.6652501833810612E-004   1.7139580609976811E-005   2.0807070657108944E-006];
base.lagrange.p2.cons = [-1.9895064180500688E-003  -2.7549796402737137E-005  -1.2664392504135868E-006  -8.2871500828307632E-008  -5.1791326782790748E-009];

base.lagrange.p3.l2 = [5.0758555997859178E-003   1.9545100193651420E-004   1.0219092647902650E-005   6.4432177483855777E-007   1.2504958182090930E-007];
base.lagrange.p3.cons = [-1.0681827257585930E-004  -4.1325299587890996E-007  -1.4657126679029631E-008  -6.9834021898529386E-010  -1.2975870378184595E-010];

base.lagrange.p4.l2 = [4.6168247478623603E-004   7.4480802362465711E-006   2.5357115815633659E-007   1.7893364406103350E-007   1.0353184674836628E-007];
base.lagrange.p4.cons = [-1.0935955654933771E-006  -3.5736583359202712E-008   1.6005608505231805E-008   1.7084406944700703E-008   1.7102998128848412E-008];

base.lagrange.p5.l2 = [8.3817996367918158E-005   1.2225175751355135E-006   1.4525462355190472E-007   1.7842747051143741E-007   1.5458890778607412E-007];
base.lagrange.p5.cons = [-1.7824191556603530E-006   4.8811011676974658E-008   4.0529472888906071E-008   4.0284754365682574E-008   4.0225528685766676E-008];

% Conservative Overset WITH Cell merging
cons.legendre.p1.l2 = [1.5293201901746848E-002   4.9564221436772145E-003   9.1538371977368063E-004   1.2837379332049549E-004   1.8030389665289336E-005];
cons.legendre.p1.cons = [-2.2724650347069608E-005  -3.4460239904249290E-006  -1.4607639921976334E-007   6.3081104562101586E-008   4.8634735169694743E-008];

cons.legendre.p2.l2 = [3.0896945737390142E-003   2.1003261538524851E-004   9.0486554730453349E-006   2.3416122157765218E-006   1.8506300192941094E-006];
cons.legendre.p2.cons = [-6.7681742588615590E-006  -5.2197857944996651E-007  -1.3933478048427927E-007  -8.2425469366542714E-008  -4.7419992106612074E-008];

cons.legendre.p3.l2 = [3.5475672734786961E-004   7.5467298241212780E-006   1.5976785972644576E-006   1.5688376183527270E-006   1.2970890030442586E-006];
cons.legendre.p3.cons = [-2.3806221037883657E-006  -5.4700524665385331E-009   7.0060799925442119E-008   4.9694170702052620E-008   2.9059694012234516E-008];

cons.legendre.p4.l2 = [3.4279885610601931E-005   7.6455828275363594E-007   9.0140337671899315E-007   8.7970560384114460E-007   7.2878103089097185E-007];
cons.legendre.p4.cons = [3.1969660585495774E-008  -2.6754110545035292E-008  -2.4652309506267844E-008  -1.7054864624721500E-008  -1.0001645414647875E-008];

cons.legendre.p5.l2 = [3.7791151921708296E-006   2.3238565705957143E-007   3.2114120738543293E-007   3.2288825917875591E-007   2.7121818198170947E-007];
cons.legendre.p5.cons = [2.0700809635898487E-008  -5.2920716420157277E-011   1.7918008743400549E-010   2.2556490703040666E-010   1.6042461109533335E-010];

cons.lagrange.p1.l2 = [0.13024082450049188        4.6774406157068807E-002   9.0300774360331515E-003   1.2840957421437196E-003   1.8035166800499229E-004];
cons.lagrange.p1.cons = [-2.7649152746223038E-005  -2.7689765223493090E-006  -3.1788921517339475E-007  -3.8449214145863664E-008  -4.7404423719754618E-009];

cons.lagrange.p2.l2 = [2.8212701936268764E-002   2.2292542159645934E-003   1.0492494398959134E-004   8.6344974133098522E-006   1.0178175446954731E-006];
cons.lagrange.p2.cons = [4.9440889759999607E-006   1.9092020459021697E-007   8.1015957831276353E-009   3.8862940643369370E-010   2.0663581956625876E-011];

cons.lagrange.p3.l2 = [3.5743158029595558E-003   9.1194666002291315E-005   4.1745818151669896E-006   3.0388115005107422E-007   1.1964245360698122E-007];
cons.lagrange.p3.cons = [-1.0783805194369833E-006  -3.6213727683254149E-008  -1.1140450228985799E-009  -3.4255043246389505E-011  -1.0600964550633307E-012];

cons.lagrange.p4.l2 = [3.7598382697954662E-004   4.0918310662661108E-006   1.7904925431443793E-007   1.7767426055219586E-007   1.0106529960313267E-007];
cons.lagrange.p4.cons = [-6.3650319792607490E-008   8.0163542470756965E-011   6.3642979775124786E-012   1.2950751582252451E-013   2.3314683517128287E-015];

cons.lagrange.p5.l2 = [4.2483599354368402E-005   5.8105125336658858E-007   1.3476290563270500E-007   1.7012798753092814E-007   1.4489730018112734E-007];
cons.lagrange.p5.cons = [2.7324858098420890E-008   1.2396872417497207E-010   6.6363581296968732E-013   4.1078251911130792E-015  -8.3266726846886741E-016];

%% Plots
close all; 

figure(10)
subplot(1,2,1);
loglog(dx,abs(base.lagrange.p1.l2),'b--',dx,cons.lagrange.p1.l2,'b-','linewidth',3)
hold on; 
loglog(dx,abs(base.lagrange.p2.l2),'r--',dx,cons.lagrange.p2.l2,'r-','linewidth',3)
loglog(dx,abs(base.lagrange.p3.l2),'g--',dx,cons.lagrange.p3.l2,'g-','linewidth',3)
loglog(dx,abs(base.lagrange.p4.l2),'m--',dx,cons.lagrange.p4.l2,'m-','linewidth',3); 
loglog(dx,abs(base.lagrange.p5.l2),'c--',dx,cons.lagrange.p5.l2,'c-','linewidth',3); 
loglog(dx,(base.lagrange.p1.l2(2)/dx(2).^2)*dx.^2,'k--')
loglog(dx,(base.lagrange.p2.l2(2)/dx(2).^3)*dx.^3,'k--')
loglog(dx,(base.lagrange.p3.l2(2)/dx(2).^4)*dx.^4,'k--')
loglog(dx,(base.lagrange.p4.l2(2)/dx(2).^5)*dx.^5,'k--')
loglog(dx,(base.lagrange.p5.l2(2)/dx(2).^6)*dx.^6,'k--')
grid on;
legend('base p=1','cons p=1','base p=2','cons p=2','base p=3','cons p=3','base p=4','cons p=4', 'base p=5','cons p=5','location','southeast');
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('L2 Error','fontsize',18)
set(gca,'fontsize',18')
ylim([1e-8 1])
title({'Lagrange','90% Cut'},'FontSize',30')

subplot(1,2,2)
loglog(dx,abs(base.lagrange.p1.cons),'b--',dx,abs(cons.lagrange.p1.cons),'b-','linewidth',3)
hold on; 
loglog(dx,abs(base.lagrange.p2.cons),'r--',dx,abs(cons.lagrange.p2.cons),'r-','linewidth',3)
loglog(dx,abs(base.lagrange.p3.cons),'g--',dx,abs(cons.lagrange.p3.cons),'g-','linewidth',3)
loglog(dx,abs(base.lagrange.p4.cons),'m--',dx,abs(cons.lagrange.p4.cons),'m-','linewidth',3); 
loglog(dx,abs(base.lagrange.p5.cons),'c--',dx,abs(cons.lagrange.p5.cons),'c-','linewidth',3); 
% loglog(dx,abs(base.lagrange.p1.cons(2)/dx(2).^2)*dx.^2,'k--')
% loglog(dx,abs(base.lagrange.p2.cons(3)/dx(3).^3)*dx.^3,'k--')
% loglog(dx,abs(base.lagrange.p3.cons(3)/dx(3).^4)*dx.^4,'k--')
% loglog(dx,abs(base.lagrange.p4.cons(2)/dx(2).^5)*dx.^5,'k--')
% loglog(dx,abs(base.lagrange.p5.cons(2)/dx(2).^6)*dx.^6,'k--')
grid on;
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('Conservation Error','fontsize',18)
ylim([1e-17 1e-2])
set(gca,'fontsize',18')
title({'Lagrange','90% Cut'},'FontSize',30')
set(gcf,'Position',scrsz);
saveas(gcf,'ConsOverset_CellMerging_Lagrange_90pctCut.png');
% close(gcf);

figure(2)
subplot(1,2,1);
loglog(dx,abs(base.legendre.p1.l2),'b--',dx,cons.legendre.p1.l2,'b','linewidth',3)
hold on; 
loglog(dx,abs(base.legendre.p2.l2),'r--',dx,cons.legendre.p2.l2,'r','linewidth',3)
loglog(dx,abs(base.legendre.p3.l2),'g--',dx,cons.legendre.p3.l2,'g','linewidth',3)
loglog(dx,abs(base.legendre.p4.l2),'m--',dx,cons.legendre.p4.l2,'m','linewidth',3); 
loglog(dx,abs(base.legendre.p5.l2),'c--',dx,cons.legendre.p5.l2,'c','linewidth',3); 
loglog(dx,(base.legendre.p1.l2(2)/dx(2).^2)*dx.^2,'k--')
loglog(dx,(base.legendre.p2.l2(2)/dx(2).^3)*dx.^3,'k--')
loglog(dx,(base.legendre.p3.l2(2)/dx(2).^4)*dx.^4,'k--')
loglog(dx,(base.legendre.p4.l2(2)/dx(2).^5)*dx.^5,'k--')
loglog(dx,(base.legendre.p5.l2(2)/dx(2).^6)*dx.^6,'k--')
grid on;
legend('base p=1','cons p=1','base p=2','cons p=2','base p=3','cons p=3','base p=4','cons p=4', 'base p=5','cons p=5','location','southeast');
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('L2 Error','fontsize',18)
set(gca,'fontsize',18')
% ylim([1e-8 1])
title({'Legendre','90% Cut'},'FontSize',30')

subplot(1,2,2)
loglog(dx,abs(base.legendre.p1.cons),'b--',dx,abs(cons.legendre.p1.cons),'b','linewidth',3)
hold on; 
loglog(dx,abs(base.legendre.p2.cons),'r--',dx,abs(cons.legendre.p2.cons),'r','linewidth',3)
loglog(dx,abs(base.legendre.p3.cons),'g--',dx,abs(cons.legendre.p3.cons),'g','linewidth',3)
loglog(dx,abs(base.legendre.p4.cons),'m--',dx,abs(cons.legendre.p4.cons),'m','linewidth',3); 
loglog(dx,abs(base.legendre.p5.cons),'c--',dx,abs(cons.legendre.p5.cons),'c','linewidth',3); 
loglog(dx,abs(base.legendre.p1.cons(2)/dx(2).^2)*dx.^2,'k--')
loglog(dx,abs(base.legendre.p2.cons(3)/dx(3).^3)*dx.^3,'k--')
loglog(dx,abs(base.legendre.p3.cons(3)/dx(3).^4)*dx.^4,'k--')
loglog(dx,abs(base.legendre.p4.cons(2)/dx(2).^5)*dx.^5,'k--')
loglog(dx,abs(base.legendre.p5.cons(2)/dx(2).^6)*dx.^6,'k--')
ylim([1e-17 1e-2])
grid on;
grid minor; 
xlabel('1/N_{elem}','fontsize',18)
ylabel('Conservation Error','fontsize',18)
set(gca,'fontsize',18')
title({'Legendre','90% Cut'},'FontSize',30')
set(gcf,'Position',scrsz);
saveas(gcf,'ConsOverset_CellMerging_Legendre_90pctCut.png');