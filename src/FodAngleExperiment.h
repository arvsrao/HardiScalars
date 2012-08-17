//Angle between Lobes of Crossing FODs is varied. The a two tensor model was //used to generate signal.    Exp(- .5x A x^t) + Exp(- .5x UBU^t x^t). U is a rotation 
// matrix
//
//    [  1           ]
// A= [       1      ]     --- Tensor for the Primary Peak.
//    [            20]
//
//       [  10           ]
//  B =  [        1      ]       ---Initial tensor for the Secondary Peak
//       [              1]
//
//
 //Primary lobe is oriented along the (0,0,1)-axis and the secondary lobe is moving

float firstPeak[10][3] = {
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
{0,0,1},
};

float secondPeak[10][3] ={
{0,1,0},
{0,0.9876883406, -0.1564344651},
{0,0.9510565163, -0.3090169944},
{0,0.8910065242, -0.4539904998},
{0,0.8090169943, -0.5877852524},
{0,0.7071067811, -0.7071067813},
{0,0.5877852522, -0.8090169944},
{0,0.4539904995, -0.8910065243},
{0,0.3090169938, -0.9510565165},
{0,0.1564344652, -0.9876883406},
};

float crossesExp[10][45] =
{

//90 degree difference between peaks. 
{ .215135112715, -.138347454587, 0.514555567889e-2, .159017032861, 0.446703547721e-2, -0.129386647117e-2, 0.559076262949e-1, -0.106967351825e-2, 0.954496495555e-1, 0.395130474129e-2, .242378639151, 0.730277468496e-2, -0.131838597047e-2, 0.649956330654e-3, 0.680452224474e-4, -0.219148302473e-1, 0.72325496736e-2, -0.111296654782e-1, -0.589045358107e-2, -0.649870714718e-1, 0.213547568132e-1, .137839280612, 0.107430551185e-1, -0.295028234242e-2, -0.114359290371e-2, -0.142582773593e-2, -0.107061470176e-2, -0.430853390621e-2, 0.249145165128e-1, -0.262779335549e-2, 0.27100300903e-1, -0.108865393957e-1, 0.167284347811e-1, -0.34165852826e-2, 0.604653485857e-1, 0.925191540491e-2, .141032920668, 0.87451570187e-2, -0.43584541305e-2, -0.393662921584e-3, 0.608545849366e-3, 0.28809306923e-2, 0.672804704073e-3, 0.578335876783e-3, -0.168364343593e-2},

//81 degree difference between peaks.
{0.216470621923,-0.152475658143,0.00969570452417,0.164377569222,0.0198912041119,-0.00208690826639,0.0600719687266,0.00130739546487,0.087202120661,0.0063913537367,0.236181260547,-0.0398662759844,-0.00183647871451,-0.028946660249,0.00080424381942,-0.0560855040653,0.00728965830418,-0.00809950032003,-0.00488226295191,-0.0689871201732,0.0212699322013,0.145537332625,0.0105138426147,-0.00689532344961,0.0146622448449,0.000365900743745,0.013040650329,-0.00796424346176,0.0386294546411,-0.0031204068524,0.0283855442299,-0.00780836355253,0.00812801157061,-0.00541780726755,0.0218526009237,0.0171915420277,0.127888084491,-0.0335141574635,-0.00944339294289,-0.0156553447707,-0.00109046844713,-0.0224575151036,0.00303046040284,-0.0263494615173,0.00245162599743},

//72 degree difference between peaks
{ .22115913682, -.185636429884, -0.526043070273e-4, .215477719049, .113249569912, -0.370362063851e-2, 0.561219413805e-1, -0.302258220944e-2, 0.352033209012e-1, -0.113948859632e-2, .203580271218, -0.597412700166e-1, -0.113928783229e-2, -0.757987026806e-1, 0.17003307657e-2, -0.631067180254e-1, 0.120255942544e-1, 0.608162475826e-1, -0.168228826107e-2, -0.533253503535e-2, -0.776812176342e-2, .225248351145, 0.569843785292e-1, 0.404937111678e-2, 0.190370535917e-1, 0.799778000795e-2, 0.841029527932e-1, 0.227184290987e-3, 0.656520050016e-1, -0.673863948313e-2, -0.206900860864e-1, 0.27107046186e-2, -0.435188788957e-1, 0.25435284682e-2, -0.33454054706e-1, 0.993056288444e-3, .114474562236, 0.266714979468e-1, 0.120447718618e-2, 0.437188798445e-1, -0.356069963127e-2, -0.324271313324e-1, -0.876870996006e-2, -0.592238894339e-1, 0.705810050133e-2},

//63 degree difference between peaks
{ .222499994085, -.175820703228, 0.457477413464e-2, .278868754495, .143587406059, -0.831133808694e-2, 0.442374236306e-1, -0.91146221019e-2, -0.907226471754e-2, -0.196132387877e-2, .177102039769, -0.488831927081e-1, -0.391016303055e-2, -0.8335386599e-1, 0.379128577392e-2, -0.438831444519e-1, 0.333474318989e-1, 0.649193221633e-1, 0.414261208376e-2, 0.145769263761e-1, -0.510284257407e-2, .224023715571, 0.994251981819e-1, 0.857668880814e-2, -0.330148569605e-2, 0.259107360291e-1, .103285739383, -0.145513236666e-1, 0.591894518479e-1, -0.212151957105e-1, -0.424681215814e-1, -0.87313202979e-2, -0.568637600297e-1, 0.962289485988e-2, -0.772478252984e-1, 0.846542245291e-2, 0.867531476164e-1, 0.770294361206e-1, -0.635520621849e-2, 0.60749573473e-1, -0.147675357236e-1, -0.396553108429e-1, -0.1032978942e-1, -0.52727779855e-1, 0.285345756809e-1},

//54 degree difference between peaks
{ 0.224562952921,-0.164032471713,-0.00488981813172,0.352494637087,0.189627481931,0.00493664467927,0.0256343953169,0.00489523082861,-0.0659143856998,0.00335519872941,0.14385723197,0.00149676540687,0.00295118369326,-0.0742885341221,0.000136409853138,0.020001497857,-0.0155254703538,0.0859495133263,0.00428280146718,0.0357651251975,0.009337568695,0.242104475173,-0.00319651670438,-0.00696218057554,-0.145781181582,-0.0228543656604,0.0641056235301,-0.00511649034321,0.00838985540676,0.018813344277,-0.0592873191532,0.00283144803117,-0.0159698289121,-0.00483018791033,0.0668632615561,0.000216700495441,0.210053329154,0.0310343052565,0.0101607792736,0.100474091689,0.0121134712096,0.0693479465045,0.0113584319345,-0.0312329573133,-0.0143380327359},

//45 degree difference between peaks
{ .223133616922, -.149595876461, 0.442872255903e-2, .394624491452, .193462060771, 0.237713016388e-2, 0.185323058325e-1, 0.108374621217e-3, -0.867690754431e-1, 0.417021468796e-2, .151677207179, 0.334746257835e-1, 0.288829717974e-3, -0.595622983748e-1, 0.33728619014e-3, 0.969877904616e-2, 0.321186719583e-2, 0.730402504345e-1, -0.119339142928e-2, -0.607735278643e-1, 0.56294570516e-2, .174123861983, 0.486021133908e-1, 0.168621590238e-2, -.117676645573, 0.454176196915e-2, 0.34159528702e-1, 0.11401630007e-2, 0.117345176248e-1, -0.12893342745e-2, -0.345077326752e-1, -0.768174051343e-2, 0.388051967146e-1, -0.247983883205e-2, .103849963569, 0.326256412423e-2, .20053663621, -0.367109090327e-1, -0.135942722284e-2, 0.295427424495e-1, -0.228741045482e-2, 0.707967774288e-1, -0.593455716547e-2, -0.149607515092e-1, 0.669930502516e-2},

//36 degree difference between peaks
{ 0.219638121695,-0.107944800009,-0.00233023276807,0.448534263104,0.171304390385,-0.00416644404092,0.010388505754,-0.00104091898498,-0.0850484432822,-0.000984096407883,0.182781649148,0.0766052357898,-0.00280139775971,-0.0390610927971,0.000860968521053,0.0063711918362,0.00818305763808,0.0417742495925,-0.000142191099352,-0.105282153218,7.86840115638e-05,0.17025781388,0.139988805758,0.00521858219771,-0.0695484126701,0.00583162024519,0.0144231397395,0.00248380614212,0.00467752974953,-0.0082573196398,-0.0107915305939,0.00595050367785,0.0573894332384,0.000260044953119,0.0689836473453,-0.00109767135574,0.14728732602,-0.0893239095416,0.0109614879874,-0.0479919322429,0.00310948977672,0.0415020000795,-0.00201033620726,-0.0101469288558,0.000417078633113},

//27 degree difference between peaks
{ .214576927175, -0.607448008827e-1, 0.13244353902e-1, .492827245995, .151318135861, -0.148389089676e-1, 0.352747151097e-2, -0.914455531013e-2, -0.62309318393e-1, 0.254368795161e-2, .237057645395, .102383510462, -0.934211571296e-2, -0.169822768103e-1, 0.291715809026e-2, 0.328980181077e-2, 0.179980900177e-2, 0.103083178677e-1, -0.257746840828e-1, -.103986093415, 0.158057663617e-2, .206514369308, .122732578638, -0.257434516661e-1, -0.454333574758e-1, 0.11402419534e-1, 0.103055110281e-1, 0.984184820707e-2, 0.845440888883e-3, -0.415161577914e-2, 0.120484031718e-2, 0.150596168058e-1, 0.247635335583e-1, -0.155322860641e-1, -0.343540933644e-1, -0.155390921387e-1, 0.605470208555e-1, -0.6263943383e-1, 0.662275702535e-2, -0.353343742469e-1, 0.191216284275e-1, 0.599446566349e-2, -0.777727690231e-2, -0.273316012232e-2, -0.228899880818e-2},

//18 degree difference between peaks
{ .21238374853, -0.321786683305e-1, 0.743095234378e-2, .531551665902, .116225384951, 0.50839344452e-2, 0.201741071846e-2, -0.109122397268e-2, -0.378108828388e-1, 0.315825323086e-2, .284625703738, 0.934340463127e-1, 0.391668197222e-2, -0.889300045607e-2, 0.813860883414e-4, -0.295871941525e-2, -0.241100634499e-2, 0.49387913758e-2, 0.882183648421e-3, -0.788586542374e-1, 0.127624744156e-1, .256595681395, .124342356434, 0.662170379436e-2, -0.234040853832e-1, 0.14661942256e-4, 0.469430260086e-2, 0.332924691422e-2, 0.258254766575e-2, 0.157192053069e-2, -0.114380531175e-2, 0.382439531496e-2, 0.152212968958e-1, 0.239810501955e-2, -0.430120507331e-1, 0.467335144183e-2, 0.569124472053e-1, -0.905948284232e-2, 0.143207168621e-1, -0.230604370048e-1, 0.179911974856e-2, 0.191803205334e-2, -0.221518850946e-2, -0.531579892561e-3, -0.125868876018e-2},

//9 degree difference between peaks
{ .210492285962, -0.860981756277e-2, -0.262997910461e-3, .548186418758, 0.598955547757e-1, 0.258126221381e-2, 0.310728724512e-3, -0.149739225201e-2, -0.109317935946e-1, -0.710106718501e-4, .322122344161, 0.551626232702e-1, 0.254978877838e-2, -0.184018986238e-2, 0.959112189245e-5, 0.116313398472e-2, 0.771211903563e-4, -0.646094491852e-4, -0.845216515638e-4, -0.24523638343e-1, 0.576216765142e-3, .318859299024, 0.844254417461e-1, 0.697158245169e-2, -0.381093685595e-2, -0.300728207173e-4, 0.352753351363e-2, 0.290299151433e-2, -0.150821891036e-2, 0.724962536062e-4, 0.105618254061e-2, 0.124713115512e-2, 0.440363360927e-2, -0.453536191635e-2, -0.198331343054e-1, 0.213809259139e-2, 0.896153298835e-1, 0.201642192663e-1, 0.797517624976e-2, -0.661726726831e-2, -0.111488016783e-2, -0.628456555561e-2, -0.21434402933e-2, 0.317485342257e-3, 0.105219428439e-2},
};