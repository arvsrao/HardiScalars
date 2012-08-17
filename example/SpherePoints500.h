#ifndef __SpherePoints500_h
#define __SpherePoints500_h

/*****************
*  500 Equidistant Sphere (Unit Sphere) points 
*  
*  Points were produced using gradient decent of
*  total columb energy of points. The inital distribution 
*  of points were drawn from the uniform distribution on S^2.
*  This scheme is implemented with MAPLE: EquidistantSpherePoints.mw  
*
******************/

double points500X[500]={0.930591712,0.375916967,-0.388609848,-0.207039687,0.187945386,-0.630915693,-0.353023507,0.801280386,-0.30777943,0.80941199,0.714682118,-0.521033786,0.223524677,-0.318277155,-0.462402293,0.743609427,-0.799671439,0.971364547,0.531185968,0.276732994,0.135280927,0.408409136,0.203002512,-0.148445628,0.418501287,0.456948761,0.0792847582,-0.314746457,0.174660916,-0.299200078,0.342382183,0.171312655,-0.890469108,-0.778694331,-0.324308596,-0.0610715431,0.252997974,-0.798728753,-0.866272445,-0.930911544,-0.262803837,-0.260857198,0.620755662,0.330880447,-0.0674259059,-0.048745118,-0.131708739,0.666486682,0.55440332,0.330556786,0.545279281,-0.498606981,0.460406586,0.305734586,-0.188541708,0.720117043,-0.573448701,-0.286393725,0.840582348,-0.0990343489,-0.0773335692,0.996802223,0.874107289,-0.316843789,0.992501436,-0.548959969,-0.154522057,0.477172623,0.586124335,0.534925904,-0.443288665,0.878313745,-0.523566745,-0.777799498,-0.699209815,-0.213444214,0.57120925,-0.926639251,0.536738412,0.318091387,-0.305573456,0.418438643,0.017218247,-0.819338625,0.909746056,0.97485829,0.964166424,0.650804816,-0.480399607,-0.539378292,0.771315165,0.167504912,-0.352309025,0.728996493,0.77790813,0.59251947,0.0784515123,0.866293806,0.216434982,0.109612933,0.775251031,0.66375479,-0.211888804,-0.201514822,-0.406225087,-0.823559851,-0.408573335,-0.28951867,0.278304522,0.41275187,0.878355252,0.963529569,0.926115977,0.691340288,-0.988880367,-0.610189634,-0.0925508372,0.0826308723,-0.505125575,-0.695427732,-0.217542971,-0.0280251295,0.717571493,-0.152309618,-0.577363309,0.649245546,0.0383024099,-0.655113752,0.468543708,-0.0867291976,-0.0789450425,0.929064398,0.641824765,0.00192083428,-0.155536123,0.0905807504,-0.735739275,0.222769858,-0.935530637,-0.599512585,0.100899595,0.282026182,-0.338313236,0.708768927,-0.990610466,0.493606769,0.863740216,-0.428957861,-0.153498272,-0.619972898,0.149161914,-0.977392373,0.976041502,0.155626335,0.584246893,-0.28201712,-0.596827754,-0.816534462,-0.461239297,-0.356547768,0.788075176,-0.691035718,-0.326762104,-0.143094184,0.329175053,-0.806592947,-0.168957437,0.476423637,0.973568696,-0.210601833,-0.904153802,0.481312344,0.659117395,0.0431693698,0.664397484,0.530095217,0.19061857,0.823888486,0.169726215,-0.640805377,-0.459887974,0.0079956197,-0.249387031,-0.823181498,0.219468358,-0.865115366,0.745361594,-0.849696513,-0.351003028,-0.102277978,0.922527161,-0.185562993,0.237500725,-0.497053211,0.728594086,-0.909144916,-0.602090767,-0.0516678932,-0.759183411,0.353690695,0.390153727,0.422674131,0.166655744,0.180946597,-0.65592721,0.795592281,0.459408183,-0.974596001,-0.846595837,-0.713822529,-0.106909107,0.609664828,-0.635205507,0.863543352,-0.0672890531,-0.296542942,0.201557561,-0.533082204,0.312666604,0.418595149,0.397033338,0.590536051,0.825541328,-0.491638649,-0.0904598259,-0.905967214,-0.692288675,0.840851193,-0.455206445,-0.590455316,0.493303099,-0.321008098,-0.165552136,0.814791658,0.120973872,-0.337178795,-0.718665238,-0.200700046,0.496987637,0.482154282,0.86966996,-0.579647378,0.516558303,-0.0452871874,-0.210331634,0.132083707,-0.473570073,-0.714911776,0.517992556,-0.592743903,0.790904764,-0.860318356,-0.578506233,-0.851185984,0.0633897546,0.459261595,-0.870059772,-0.0493292705,-0.781698708,-0.171246846,-0.754415039,0.371525013,0.362637577,-0.708788805,-0.694256975,0.148928631,0.239905024,-0.726833627,-0.899937175,-0.052235426,0.838193696,0.552479313,-0.49562815,0.619387696,0.45568727,0.384583834,-0.238650146,-0.815207365,0.47536864,0.450184757,0.0238065595,-0.0953504811,-0.844338072,0.933078334,-0.840521101,0.309707751,-0.161288773,-0.62258257,0.1136746,0.889668109,-0.0245830199,0.77881876,0.121143709,0.305279975,-0.0603421885,-0.334926942,0.0820468126,-0.90012875,-0.779159959,-0.583228373,0.46702244,-0.519876168,-0.370826517,-0.359188183,-0.340910099,-0.451501729,-0.872081776,0.652407062,0.228235166,0.602939636,0.313817624,0.349013581,0.347472048,0.597513205,0.594486305,-0.0298138608,0.975379717,-0.79167037,-0.566519072,-0.525167328,0.547033911,-0.110836383,-0.439576176,0.626772968,-0.40846836,0.623648075,-0.929180249,0.213857329,-0.652937494,-0.394806662,-0.655111294,-0.91483469,-0.246847657,0.663713531,0.899815083,-0.00732920804,-0.961424913,-0.390713858,0.514345099,-0.976658234,-0.945173282,0.0688748336,-0.925270933,0.73191988,-0.513838998,-0.994852369,-0.382097716,0.250396733,0.448936901,0.404199177,0.684407324,-0.0449178796,0.711252418,-0.246172702,0.588534409,0.841307214,0.00990727092,-0.676813835,0.836382828,-0.463558397,0.752356314,-0.968560112,0.941492306,-0.719456278,0.000688316497,-0.95431339,-0.0721758189,-0.573309976,-0.887502671,-0.824626831,0.0871153486,-0.786543801,0.93517677,-0.418703084,0.907002731,0.692829618,0.770899609,0.77204379,0.0453112065,0.184988782,0.332590372,0.55006086,0.656390874,0.67364394,-0.458541743,-0.291644876,-0.66914204,-0.604159349,-0.62213925,-0.876812021,0.585132045,0.0571146766,-0.15375428,-0.403425382,0.340077956,0.311505536,0.882536045,0.858731538,0.520633183,0.462356659,0.396171742,-0.954560567,0.812028478,-0.459525182,0.0363293201,-0.388112692,0.204158229,0.369758929,0.801583216,0.0212021096,-0.454388356,-0.755430726,0.945500769,0.0947933296,-0.45115675,0.03676326,0.0508258417,0.708828262,-0.783448885,0.43594702,0.618837009,0.602190495,0.790051916,0.262141978,0.856733824,-0.872719633,0.303369652,0.349609539,-0.27413648,0.0511547123,0.927606237,0.0013649508,0.728239961,-0.110186068,-0.720707068,-0.485340348,0.289914232,-0.0918338178,-0.0578090452,-0.202447524,-0.628951408,-0.579628132,0.73805732,0.869899313,-0.964499389,-0.00139047538,-0.394622384,-0.242709099,0.37456313,-0.595740916,-0.348628323,0.170107424,-0.46752972,-0.249267775,-0.686113378,0.44766907,0.257295108,-0.720554192,-0.152752739,-0.747116618,0.704918891,-0.220512584,-0.18957326,0.953867209,0.102141138,-0.827457211,0.334571467,0.239199908,-0.157398532,0.913365082,0.0407373696,-0.964201776,-0.493997012,-0.712506676,-0.579703915,0.223481194,-0.447473046,0.695212952,-0.928003853,-0.288513733,0.192393264,-0.782173088,0.582133533,-0.726129067,-0.689173273,0.263363193,0.981662117,-0.742512792,0.922738237,0.0600330268,-0.927398416,0.761064735,0.0958656662,0.56550275,-0.572069854,-0.247981338,-0.997206046,0.207219069,0.996296818,0.143851098};

double points500Y[500]={0.356394467,-0.897033827,0.0202400413,-0.919527793,0.660719636,-0.112321752,-0.932888997,-0.207527393,0.224039302,0.519893723,0.220506738,-0.567351824,-0.949503823,-0.761160226,0.573605212,-0.338915396,-0.262512975,-0.23736027,0.487609983,0.68580142,-0.783056321,-0.280536694,-0.866084651,0.609246537,0.905020407,-0.642365036,-0.996323458,0.933948457,0.976396219,0.540472253,-0.433200629,-0.672333757,-0.422661265,-0.412315412,0.667178618,-0.472155081,-0.962702696,0.373492967,0.291666244,0.26954814,0.964848524,-0.538726845,-0.446161532,-0.241327851,-0.898074578,0.958317554,-0.645461287,0.741392829,0.619926603,0.911721986,0.838173229,-0.31378171,0.371172362,-0.710953518,-0.849443605,-0.475828141,0.144723282,-0.472705644,-0.499320466,-0.00316031447,0.824537266,-0.0202271008,0.289507189,0.763372623,-0.0986279162,0.771939164,0.473970175,0.422865545,-0.181950831,0.783315495,0.58072148,-0.474495069,-0.0187889613,-0.570549342,-0.35731274,-0.96354032,0.269097079,-0.284489702,0.16065328,-0.573380179,0.398686931,-0.752714854,0.991161172,-0.45386419,-0.244997579,0.222672962,0.205602715,0.0923646962,-0.808730342,-0.454545289,0.597900022,-0.660702788,-0.31824852,0.442827664,0.41334093,0.742782891,-0.945021524,0.499212312,0.911762402,0.993721018,0.628588493,-0.737904138,-0.405211804,-0.963667569,-0.52713966,-0.00784016779,0.791466934,-0.662066828,0.172814862,0.687070506,0.0877970686,-0.114812015,-0.263010153,0.718287405,0.142389696,-0.739527314,0.898910542,-0.926887135,-0.225919105,0.374183056,0.929734853,-0.827876583,-0.34546424,-0.569624282,0.451005676,-0.291302793,-0.0295258282,0.00946113269,0.035137505,0.996222212,0.0707103049,-0.34490989,0.381770726,0.557141502,-0.0721841007,0.573719963,-0.677202292,-0.936735043,-0.133864382,-0.665791462,-0.264254819,0.932416787,-0.303562525,0.659488375,-0.136709537,-0.859972942,-0.389184008,0.117867258,0.359986544,-0.776999856,-0.138385238,-0.132295875,-0.160809881,0.836723372,0.476478793,0.775096663,0.66336325,0.576501336,0.689795079,-0.929743873,0.605905226,0.436158609,0.931480578,0.730827795,0.419250557,-0.546175126,0.97547877,0.748316223,0.0236278671,0.932301828,-0.427190699,-0.438687405,0.104413389,0.0776560599,0.706879636,-0.833120609,-0.778973747,-0.385993616,0.343559574,0.0932346566,0.269632794,0.370216392,-0.0647285181,-0.559095309,0.757225162,-0.136662248,-0.00221498057,0.100750657,-0.900501104,0.169149,0.306133465,0.620568076,-0.897852944,0.848949989,0.556113044,0.0076076986,0.787687891,-0.357490438,0.059024171,-0.789202477,-0.837574796,0.875134793,0.283856556,-0.528124758,-0.671479444,0.131254996,-0.299584721,-0.145820965,-0.201269039,-0.612903058,-0.69267575,-0.741161543,0.752014857,0.47274247,-0.995338509,0.6663026,-0.546591998,0.765942664,-0.695955489,0.902665698,0.509464866,-0.583742971,0.402997452,-0.860835151,-0.803765798,0.31224469,0.560882343,-0.0565787246,0.421064849,0.306491361,0.630523483,0.372762401,0.821313832,-0.577836216,0.965379798,0.535756609,0.594290779,-0.920347062,0.741225657,-0.513984229,0.472534808,0.251835611,0.13329268,-0.313484425,-0.977610918,0.74654152,-0.831598871,-0.127485941,-0.737610062,-0.159382704,-0.519597392,0.419967345,0.672650663,0.146701522,0.916204037,-0.133450485,0.247314034,0.958961658,0.238714731,-0.25180364,0.640639535,0.844590169,-0.438457387,-0.523722954,0.719504236,0.622152925,0.96231138,-0.0143137628,0.4355104,-0.960688652,0.438252486,-0.708117136,-0.868166195,0.54182051,-0.685267349,0.83508759,-0.223865621,-0.0709591179,0.841535584,-0.105497824,-0.601007921,-0.514739222,0.520460474,0.348544402,-0.435148471,0.557739434,-0.757011764,-0.271449088,0.952793263,-0.377650705,-0.179207678,-0.618753346,0.112200995,0.790981717,-0.992786302,-0.844858512,-0.430188373,-0.0322605575,0.554269481,-0.568445324,-0.827512038,-0.74908477,-0.0259260191,0.875959154,-0.910808954,0.312707357,-0.308926603,-0.694046465,-0.359766456,-0.328402396,-0.125549858,-0.934227626,0.119327012,-0.0262376415,0.30312045,0.658093147,0.162651698,0.571032972,0.82352004,-0.0446311644,-0.211122896,0.898721582,0.784231855,-0.761500159,0.911321441,0.667951553,-0.151391518,-0.411440981,0.707949355,-0.762430772,-0.506647821,0.380773358,0.100314171,-0.118572941,0.349286938,0.815966147,0.272539905,-0.146421138,-0.372629644,0.153016254,0.301008672,-0.985499238,0.11801257,-0.632400667,0.119820371,0.00276611826,0.190224617,0.838778587,0.282635099,0.240643579,-0.596341014,0.512355462,-0.701797212,-0.747974294,-0.599377107,-0.233528676,0.219443052,0.665832286,-0.0546078332,-0.385916355,-0.628090083,0.0156036648,0.0353973519,-0.553528807,0.69223494,-0.285464951,-0.966664195,-0.429961081,-0.428197135,0.454950892,0.842067221,0.437203782,0.153850431,-0.557858799,-0.101289831,0.374235593,0.300494404,-0.516233779,0.902506142,0.171720268,-0.869981548,0.820893097,-0.636450236,0.527235616,0.682986713,-0.154778893,-0.426899416,-0.763804815,-0.780839975,0.446770654,0.795496638,-0.847891266,0.213905782,-0.665993702,-0.583722888,0.270981907,-0.201728298,0.253025963,-0.801054631,-0.571960629,0.0506007355,-0.288774743,-0.355157064,0.883719677,-0.60273792,-0.435576271,-0.977046381,-0.908494499,0.262527119,0.421283279,-0.655402232,0.144815315,-0.112158212,-0.977049274,-0.171744538,0.252428224,0.750389839,-0.485791295,0.297885199,-0.898841715,0.646967786,-0.471570864,0.0825175678,-0.272018603,-0.488734573,-0.286971466,0.371378999,0.639525413,-0.62884901,-0.0970448575,0.195844505,-0.724640929,-0.0411445574,0.333763674,-0.671307619,0.842495284,0.956836053,-0.167300501,-0.910485823,0.490612144,-0.302915559,0.415874413,0.564562363,0.13061049,-0.00313529218,0.991718985,0.868606968,0.859671105,0.743239369,0.564652439,-0.850704775,0.488023846,-0.742609602,0.0684282761,0.284716285,0.564450299,0.527213887,-0.674684676,0.737368018,-0.262357464,-0.192035576,-0.843216703,-0.373153461,-0.254608242,-0.876131106,-0.557611884,-0.0528648704,0.0109892996,0.976667949,-0.00092614433,-0.738310634,0.162127928,-0.848989577,0.501498994,0.556218657,0.894054227,0.453804702,0.223283506,0.144967764,0.857201963,0.0144446109,-0.374351842,-0.813088454,0.676235231,0.194466897,-0.818439189,0.0393719608,-0.167540524,-0.380661658,-0.46338146,-0.290187418,-0.179526511,-0.29342473,-0.0224112063,-0.643761442,0.269501237,0.0130792024,-0.154078443,0.0781968643,0.436218842};

double points500Z[500] = {-0.0835586586,0.232415033,-0.921180073,0.33407066,0.726722846,-0.767677803,-0.0713619161,-0.561143586,0.924704393,0.273061802,-0.663781778,-0.637679937,0.220179947,0.565097127,0.676136954,0.57635178,-0.5400116,-0.0105365457,0.692876592,-0.673123512,0.607060021,0.868620252,-0.456823113,-0.778962485,-0.0761231566,0.615276353,-0.0324575788,-0.169336788,-0.127058994,-0.786364455,-0.833733564,0.720152271,0.168588916,-0.472896542,0.670594234,-0.879397433,-0.0958934039,0.47173656,-0.405589514,-0.246470073,-0.00121213484,0.801078591,0.644672238,-0.912293263,-0.434644452,-0.281516214,0.752351338,-0.0783082105,-0.555272875,0.243916444,0.0116680908,0.80804215,0.806385053,0.633302027,0.492846413,0.504994083,0.806357091,-0.83338347,-0.210000922,-0.995078997,0.560497829,-0.0773058359,-0.390028249,-0.562914072,0.072204114,-0.3205509,-0.866876696,-0.770383683,-0.789526541,-0.316656142,-0.682830668,-0.0584755929,0.851777458,-0.263631162,-0.619219864,0.161343169,0.775439717,-0.245774912,-0.828312984,-0.75501857,-0.86468121,0.508261203,-0.131541112,0.350273485,0.335169061,-0.00825022458,0.167662252,0.753605902,0.339369197,0.708843874,-0.21814784,-0.731719844,-0.880111487,0.521984457,-0.473295062,0.311759929,0.317458153,-0.0179473851,0.349063349,0.0224397551,-0.0621477717,-0.122176359,0.889329258,-0.17531855,-0.746394639,-0.567175197,0.454585441,0.691264388,0.944818298,-0.597971583,0.469876288,0.241720853,-0.270434569,0.078178074,0.0429044809,0.28419705,0.428245467,-0.36613164,-0.832951806,0.613487824,-0.297099577,0.560209743,-0.604769056,-0.80766946,0.6806214,0.702583023,0.998829891,0.755471018,-0.88274127,-0.00439891467,0.994368007,-0.133703078,0.665065473,-0.830415352,0.985189307,0.814027353,0.0092073613,-0.270001942,-0.326898692,0.444191839,-0.959160394,-0.225964969,0.89072664,0.250442992,-0.00118667781,0.129610556,-0.320138482,0.895601732,0.920243429,0.109109252,0.979081329,-0.164930141,-0.14657138,-0.525047087,-0.656980615,0.565412687,0.451382134,0.0302932458,0.558069046,0.0919239941,0.108703146,-0.576398563,0.159907661,-0.667394028,-0.84609264,0.226054307,0.141048058,-0.461566191,-0.227168918,0.294041101,-0.00374289206,-0.758875344,-0.744756406,-0.99604515,-0.242687379,-0.157826205,-0.597381337,0.414989997,-0.923666515,-0.762020845,0.846050358,0.928911133,-0.966238236,-0.098917428,0.615177774,-0.482595932,0.666656724,0.517556896,0.256699507,-0.980269265,-0.235001996,0.761880331,0.370747213,-0.179504936,0.399861151,-0.416410428,0.130516264,0.932486469,0.648194951,-0.502058108,0.382424542,-0.235553972,0.944272905,-0.82966413,-0.344782324,0.591443191,-0.836177683,0.169996709,0.492713165,0.33883217,-0.713281676,-0.281048687,-0.176033004,-0.175520473,0.0690900594,-0.684180626,0.812780252,0.359381691,-0.64644068,0.0998636032,0.763419989,0.557235423,0.395062741,-0.131355555,-0.588028539,0.285878751,-0.454034567,-0.538301141,0.784532654,0.74660938,0.599243014,0.870633099,-0.545926808,0.0471154012,0.23109991,0.774128747,0.361024578,-0.335679574,0.451207062,0.709469845,0.142776089,-0.774975962,0.845813562,-0.948512828,-0.00612344225,-0.652096352,-0.290128771,0.687494354,-0.433145596,0.789462994,0.323246352,-0.288928635,-0.461380129,-0.503945516,0.395666402,-0.878219651,0.426417358,0.279211679,0.576161789,0.952506926,-0.142964802,-0.385534189,0.822343668,0.47258089,-0.0179138847,-0.768600027,0.128071808,-0.686664398,0.0210659231,-0.272669717,-0.324601426,0.439700728,0.0253020748,0.568145595,-0.568117745,0.393349707,-0.944949888,0.574805842,0.256597966,0.886681393,-0.798888432,0.852028062,-0.127334654,0.0887785038,-0.322753912,-0.770070018,0.633181744,-0.733966169,-0.281537001,0.256652684,0.983504084,-0.102886524,0.986273359,-0.53023774,-0.103605885,0.417178663,-0.899002939,0.434427772,-0.292737255,0.580271124,-0.31163098,-0.410610248,0.928340205,0.321992871,-0.232824729,-0.835679542,0.379523029,0.304408492,-0.90469708,-0.7270596,0.941145967,0.0735408916,-0.930066793,-0.801429695,-0.74478455,0.752346026,-0.148925596,0.217208125,-0.0295107705,-0.849827828,0.810050012,-0.424281526,-0.437896316,0.165145857,-0.0514473499,-0.40608349,0.337201236,0.885991739,0.26922173,-0.512666576,-0.560479418,-0.134495871,-0.963848174,0.738528812,0.261402851,-0.578053224,-0.0372040004,-0.908792348,0.772396445,-0.150747206,0.126654038,0.155072593,0.360481254,-0.253698022,-0.849477818,0.101297148,-0.904331759,0.483479015,-0.847686888,0.882447558,0.419480643,0.857597962,0.0400083912,-0.616387424,-0.542562745,0.487510543,0.975575006,-0.313991082,0.545418875,-0.797610293,0.198652523,-0.248289822,0.3351705,-0.419509861,-0.721671888,-0.0882933476,0.245664376,-0.69745906,-0.170253407,0.336170009,0.532291002,-0.436121201,-0.319021088,0.716576086,0.408762054,-0.616388547,0.561619894,-0.370743943,-0.428286772,-0.967621465,-0.36402713,-0.15351734,-0.405070509,-0.517905828,-0.568567164,0.94392095,0.608281035,-0.227142436,-0.0568479373,0.177754501,0.157497903,-0.527084731,0.964678152,-0.627455456,0.737302228,-0.910787081,-0.424777381,0.445598482,0.295385454,-0.677559856,0.916781063,0.0736428634,-0.463112524,0.0887476154,0.797111774,0.812183385,0.0608258704,-0.194720514,-0.537162786,-0.906681239,0.603306834,-0.639024994,-0.305694015,-0.190758854,0.875763324,-0.966916984,0.65903856,-0.511438475,-0.54541017,-0.0450973126,0.44550358,-0.644195256,-0.607460963,0.925898193,0.16475944,-0.394978253,0.877521791,0.684675409,-0.727597494,-0.993964532,0.31810627,0.689125286,-0.684086021,-0.936194873,0.172995386,0.233765789,-0.0203593999,-0.981619627,0.409479768,0.847534498,0.716004393,-0.700770791,-0.369514182,-0.475621788,0.264066465,0.128419322,-0.299658321,0.449508492,0.55434439,0.571192073,-0.393395067,-0.856093564,0.47952783,0.966014051,-0.669466257,-0.693533174,0.80984242,-0.160007015,0.657993166,0.610725241,-0.682796964,-0.490265022,-0.908195174,0.159097425,0.471127874,0.0662083831,-0.940886411,0.970908152,-0.146100715,-0.407140478,-0.673229436,0.209832005,0.187573056,0.490747384,-0.595452917,-0.388230852,-0.770603118,0.683244778,-0.343210134,-0.426573113,-0.98121159,0.498062203,0.00277757488,0.124267815,0.698013485,0.510682996,0.186518997,-0.648540612,0.0604222448,0.884123101,0.236058129,0.623338349,0.951163341,0.824441858,-0.5082394,-0.930523691,-0.0735461437,-0.966084929,0.0357477239,0.888267968};

#endif