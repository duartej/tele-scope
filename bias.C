
// Daniel Pitzl, Dec 2018
// bias scan: eff

// root -l bias.C

//void bias()
{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.22 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(111111);
  gStyle->SetStatFormat("8.6g"); // more digits, default is 6.4g
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize(1); // no 'shadow'

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.90);

  gStyle->SetPalette(1); // rainbow colors

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  gROOT->ForceStyle();

  //--------------------------------------------------------------------
  // square canvas:
  //            topleft x, y, width x, y
  TCanvas c1("c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "bias.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  string HASH {"#"};

  //----------------------------------------------------------------------------
  // chip 509 = FTH150P RD53A 100x25 default  5E15 n_eq at CERN PS

  cout << "try to open dat file for c509";
  ifstream Dstream09( "bias509.dat" ); // turn 0 = vertical
  if( !Dstream09 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb09;
  vector <double> q009;
  vector <double> ef09;
  vector <double> ls09;

  while( Dstream09.good() && ! Dstream09.eof() ) {

    string rl;
    getline( Dstream09, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb09.push_back(vb);
    ef09.push_back(eff);
    ls09.push_back(100-eff);
    q009.push_back(q0);

  } // while lines

  int n09 = vb09.size();
  cout << "    " << n09 << " runs for chip 509i" << endl;

  //----------------------------------------------------------------------------
  // chip 509 = FTH150P RD53A 100x25 default  5E15 n_eq at CERN PS with turn

  cout << "try to open dat file for c509t34";
  ifstream Dstream0934( "bias509t34.dat" ); // turn 34
  if( !Dstream0934 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb0934;
  vector <double> ef0934;
  vector <double> ls0934;

  while( Dstream0934.good() && ! Dstream0934.eof() ) {

    string rl;
    getline( Dstream0934, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb0934.push_back(vb);
    ef0934.push_back(eff);
    ls0934.push_back(100-eff);

  } // while lines

  int n0934 = vb0934.size();
  cout << "    " << n0934 << " runs for chip 509i" << endl;

  //----------------------------------------------------------------------------
  // chip 509 = FTH150P RD53A 100x25 default  5E15 n_eq at CERN PS
  // annealed: Jun 2020

  cout << "try to open dat file for c509";
  ifstream Dstream09a( "bias509n.dat" ); // turn 0 = vertical
  if( !Dstream09a ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb09a;
  vector <double> q009a;
  vector <double> ef09a;
  vector <double> ls09a;

  while( Dstream09a.good() && ! Dstream09a.eof() ) {

    string rl;
    getline( Dstream09a, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb09a.push_back(vb-0.02*ib); // 2020: 0.02 MOhm drop
    ef09a.push_back(eff);
    ls09a.push_back(100-eff);
    q009a.push_back(q0);

  } // while lines

  int n09a = vb09a.size();
  cout << "    " << n09a << " runs for chip 509i" << endl;

  //----------------------------------------------------------------------------
  // chip 512 = FTH150P RD53A 100x25 bias dot  5.6E15 n_eq at CERN PS

  cout << "try to open dat file for c512";
  ifstream Dstream12( "bias512.dat" ); // turn 0 = vertical
  if( !Dstream12 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb12;
  vector <double> ef12;
  vector <double> ls12;

  while( Dstream12.good() && ! Dstream12.eof() ) {

    string rl;
    getline( Dstream12, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb12.push_back(vb);
    ef12.push_back(eff);
    ls12.push_back(100-eff);

  } // while lines

  int n12 = vb12.size();
  cout << "    " << n12 << " runs for chip 512i" << endl;

  //----------------------------------------------------------------------------
  // chip 512 = FTH150P RD53A 100x25 bias dot  5.6E15 n_eq at CERN PS turn

  cout << "try to open dat file for c512turn";
  ifstream Dstream1234( "bias512turn.dat" ); // turn 34
  if( !Dstream1234 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb1234;
  vector <double> ef1234;
  vector <double> ls1234;

  while( Dstream1234.good() && ! Dstream1234.eof() ) {

    string rl;
    getline( Dstream1234, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb1234.push_back(vb);
    ef1234.push_back(eff);
    ls1234.push_back(100-eff);

  } // while lines

  int n1234 = vb1234.size();
  cout << "    " << n1234 << " runs for chip 512i" << endl;

  //----------------------------------------------------------------------------
  // chip 511 = FDB150P RD53A 50x50 bias dot with 7E15 CERN PS p/cm2 on Lin

  cout << "try to open dat file for c511turn";
  ifstream Dstream11t( "bias511turn.dat" ); // turn 34
  if( !Dstream11t ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb11t;
  vector <double> ef11t;
  vector <double> ls11t;

  while( Dstream11t.good() && ! Dstream11t.eof() ) {

    string rl;
    getline( Dstream11t, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb11t.push_back(vb);
    ef11t.push_back(eff);
    ls11t.push_back(100-eff);

  } // while lines

  int n11t = vb11t.size();
  cout << "    " << n11t << " runs for chip 511i" << endl;

  //----------------------------------------------------------------------------
  // chip 521 = FTH150P RD53A 100x25 bias dot wiggle  5.3E15 n_eq at KaZ

  cout << "try to open dat file for c521";
  ifstream Dstream21( "bias521.dat" ); // turn 18
  if( !Dstream21 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb21;
  vector <double> ef21;
  vector <double> ls21;
  vector <double> q021;

  while( Dstream21.good() && ! Dstream21.eof() ) {

    string rl;
    getline( Dstream21, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb21.push_back(vb);
    ef21.push_back(eff);
    ls21.push_back(100-eff); // loss = ineff
    q021.push_back(q0);

  } // while lines

  int n21 = vb21.size();
  cout << "    " << n21 << " runs for chip 521i" << endl;

  //----------------------------------------------------------------------------
  // chip 521 = FTH150P RD53A 100x25 bias dot wiggle  5.3E15 n_eq at KaZ
  // annealed

  cout << "try to open dat file for c521";
  ifstream Dstream21a( "bias521a.dat" ); // turn 18
  if( !Dstream21a ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb21a;
  vector <double> ef21a;
  vector <double> ls21a;
  vector <double> q021a;

  while( Dstream21a.good() && ! Dstream21a.eof() ) {

    string rl;
    getline( Dstream21a, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb21a.push_back(vb-0.2*ib); // Jun 2020: 0.02 MOhm
    ef21a.push_back(eff);
    ls21a.push_back(100-eff); // loss = ineff
    q021a.push_back(q0);

  } // while lines

  int n21a = vb21a.size();
  cout << "    " << n21a << " runs for chip 521ai" << endl;

  //----------------------------------------------------------------------------
  // chip 542 = FDB150P RD53A100x25 P1 default KAZ 7.4E15 neq/cm^2

  cout << "try to open dat file for c542";
  ifstream Dstream42( "bias542.dat" );
  if( !Dstream42 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb42;
  vector <double> ef42;
  vector <double> ls42;
  vector <double> q042;

  while( Dstream42.good() && ! Dstream42.eof() ) {

    string rl;
    getline( Dstream42, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb42.push_back(vb);
    ef42.push_back(eff);
    ls42.push_back(100-eff); // ineff
    q042.push_back(q0);

  } // while lines

  int n42 = vb42.size();
  cout << "    " << n42 << " runs for chip 542i" << endl;

  //----------------------------------------------------------------------------
  // chip 542 = FDB150P RD53A100x25 P1 default KAZ 7.4E15 neq/cm^2
  // annealed

  cout << "try to open dat file for c542";
  ifstream Dstream42a( "bias542a.dat" );
  if( !Dstream42a ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb42a;
  vector <double> ef42a;
  vector <double> ls42a;
  vector <double> q042a;

  while( Dstream42a.good() && ! Dstream42a.eof() ) {

    string rl;
    getline( Dstream42a, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb42a.push_back(vb-0.2*ib); // 0.2 MOhm drop
    ef42a.push_back(eff);
    ls42a.push_back(100-eff); // ineff
    q042a.push_back(q0);

  } // while lines

  int n42a = vb42a.size();
  cout << "    " << n42a << " runs for chip 542ai" << endl;

  //----------------------------------------------------------------------------
  // chip 547 FDB150P RD53A 50x50 P2 = open p-stop KaZ 4.4 n_eq p/cm2

  cout << "try to open dat file for c547";
  ifstream Dstream47( "bias547.dat" );
  if( !Dstream47 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb47;
  vector <double> ef47;
  vector <double> ls47;
  vector <double> q047;

  while( Dstream47.good() && ! Dstream47.eof() ) {

    string rl;
    getline( Dstream47, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb47.push_back(vb);
    ef47.push_back(eff);
    ls47.push_back(100-eff); // ineff
    q047.push_back(q0);

  } // while lines

  int n47 = vb47.size();
  cout << "    " << n47 << " runs for chip 547i" << endl;

  //----------------------------------------------------------------------------
  // chip 563 FTH150P RD53A 100x25 P1 default design, Birmingham 1E16 n_eq/cm2
  // thr 900

  cout << "try to open dat file for c563";
  ifstream Dstream63( "bias563_0900.dat" );
  if( !Dstream63 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb63;
  vector <double> ef63;
  vector <double> ls63;
  vector <double> q063;

  while( Dstream63.good() && ! Dstream63.eof() ) {

    string rl;
    getline( Dstream63, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb63.push_back(vb-0.2*ib); // 0.2 MOhm drop
    ef63.push_back(eff);
    ls63.push_back(100-eff); // ineff
    q063.push_back(q0);

  } // while lines

  int n63 = vb63.size();
  cout << "    " << n63 << " runs for chip 563i" << endl;

  //----------------------------------------------------------------------------
  // chip 564 FTH150P RD53A 50x50 P1 default design, Birmingham 5E15 n_eq/cm2
  // thr 1000

  cout << "try to open dat file for c564";
  ifstream Dstream64( "bias564.dat" );
  if( !Dstream64 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb64;
  vector <double> ef64;
  vector <double> ls64;
  vector <double> q064;

  while( Dstream64.good() && ! Dstream64.eof() ) {

    string rl;
    getline( Dstream64, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb64.push_back(vb-0.2*ib); // 0.2 MOhm drop
    ef64.push_back(eff);
    ls64.push_back(100-eff); // ineff
    q064.push_back(q0);

  } // while lines

  int n64 = vb64.size();
  cout << "    " << n64 << " runs for chip 564i" << endl;

  //----------------------------------------------------------------------------
  // chip 565 FTH150P RD53A 100x25 P1 default design, Birmingham 5E15 n_eq/cm2
  // thr 940

  cout << "try to open dat file for c565";
  ifstream Dstream65( "bias565.dat" );
  if( !Dstream65 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb65;
  vector <double> ef65;
  vector <double> ls65;
  vector <double> q065;

  while( Dstream65.good() && ! Dstream65.eof() ) {

    string rl;
    getline( Dstream65, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb65.push_back(vb-0.2*ib); // 0.2 MOhm drop
    ef65.push_back(eff);
    ls65.push_back(100-eff); // ineff
    q065.push_back(q0);

  } // while lines

  int n65 = vb65.size();
  cout << "    " << n65 << " runs for chip 565i" << endl;

  //----------------------------------------------------------------------------
  // chip 577 FDB150P RD53A 100x25 P1 default design, KaZ 1E16 n_eq/cm2 on Zh card
  // thr 1160

  cout << "try to open dat file for c577";
  ifstream Dstream77( "bias577.dat" );
  if( !Dstream77 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb77;
  vector <double> ef77;
  vector <double> ls77;
  vector <double> q077;

  while( Dstream77.good() && ! Dstream77.eof() ) {

    string rl;
    getline( Dstream77, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb77.push_back(vb-0.02*ib); // 0.02 MOhm drop
    ef77.push_back(eff);
    ls77.push_back(100-eff); // ineff
    q077.push_back(q0);

  } // while lines

  int n77 = vb77.size();
  cout << "    " << n77 << " runs for chip 577i" << endl;

  //----------------------------------------------------------------------------
  // chip 578 FTH150P RD53A 100x25 P1 default design, KaZ 5E15 n_eq/cm2
  // thr 1300

  cout << "try to open dat file for c578";
  ifstream Dstream78( "bias578.dat" );
  if( !Dstream78 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb78;
  vector <double> ef78;
  vector <double> ls78;
  vector <double> q078;

  while( Dstream78.good() && ! Dstream78.eof() ) {

    string rl;
    getline( Dstream78, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb78.push_back(vb-0.02*ib); // 0.02 MOhm drop
    ef78.push_back(eff);
    ls78.push_back(100-eff); // ineff
    q078.push_back(q0);

  } // while lines

  int n78 = vb78.size();
  cout << "    " << n78 << " runs for chip 578i" << endl;

  //----------------------------------------------------------------------------
  // chip 589 FTH150P RD53A 100x25 P1 default design, KaZ 2E16 n_eq/cm2
  // thr 1170

  cout << "try to open dat file for c589";
  ifstream Dstream89( "bias589.dat" );
  if( !Dstream89 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb89;
  vector <double> ef89;
  vector <double> ls89;
  vector <double> q089;

  while( Dstream89.good() && ! Dstream89.eof() ) {

    string rl;
    getline( Dstream89, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;

    vb89.push_back(vb-0.02*ib); // 0.02 MOhm drop
    ef89.push_back(eff);
    ls89.push_back(100-eff); // ineff
    q089.push_back(q0);

  } // while lines

  int n89 = vb89.size();
  cout << "    " << n89 << " runs for chip 589i" << endl;

  //----------------------------------------------------------------------------
  // chip 589 FTH150P RD53A 100x25 P1 default design, KaZ 2E16 n_eq/cm2
  // thr 1240, LDAC_LIN: 200, KRUM_CURR_LIN: 20

  cout << "try to open dat file for c589b";
  ifstream Dstream89b( "bias589b.dat" );
  if( !Dstream89b ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb89b;
  vector <double> ef89b;
  vector <double> ls89b;
  vector <double> q089b;
  vector <double> q089r;
  vector <double> msk89;

  while( Dstream89b.good() && ! Dstream89b.eof() ) {

    string rl;
    getline( Dstream89b, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;
    int run;
    double eff;
    double q0;
    int msk;

    strm >> vb;
    strm >> ib;
    strm >> run;
    strm >> eff;
    strm >> q0;
    strm >> msk;

    vb89b.push_back(vb-0.02*ib); // 0.02 MOhm drop
    ef89b.push_back(eff);
    ls89b.push_back(100-eff); // ineff
    q089b.push_back(q0);
    q089r.push_back(q0*0.79); // Krum_Curr_Lin 20->29
    msk89.push_back(msk);

  } // while lines

  int n89b = vb89b.size();
  cout << "    " << n89b << " runs for chip 589bi" << endl;

  //----------------------------------------------------------------------------
  // eff vs bias

  he09 = new
    TH1F( "he09",
	  "RD53A 5.2#upoint10^{15 }n_{PS}/cm^{2};sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  he09->SetStats(kFALSE); // no statistics
  he09->SetMinimum(80);
  he09->SetMaximum(100);
  he09->Draw();

  ge09 = new TGraph( n09, &vb09[0], &ef09[0] );
  ge09->SetLineColor(2);
  ge09->SetMarkerColor(2);
  ge09->SetMarkerStyle(20);
  ge09->SetMarkerSize(1.4);
  ge09->Draw("Pc"); // without axis option: overlay

  ge0934 = new TGraph( n0934, &vb0934[0], &ef0934[0] );
  ge0934->SetLineColor(4);
  ge0934->SetMarkerColor(4);
  ge0934->SetMarkerStyle(24);
  ge0934->SetMarkerSize(1.6);
  ge0934->Draw("Pc"); // without axis option: overlay

  TLegend * lgnd09 = new TLegend( 0.3, 0.2, 0.94, 0.36 );
  lgnd09->AddEntry( ge09,   "509 Lin default, turn 0 ", "p" );
  lgnd09->AddEntry( ge0934, "509 Lin default, turn 34^{o} ", "p" );
  lgnd09->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // log 1-eff vs bias

  hi = new
    TH1F( "hi",
	  "RD53A PS vs KaZ;sensor bias [V];RD53A Lin inefficiency [%]",
	  82, 0, 820 ); // axis range
  hi->SetStats(kFALSE); // no statistics
  hi->SetMinimum(1e-2);
  hi->SetMaximum(100);
  c1.SetGrid();
  gPad->SetLogy(1);
  hi->Draw();

  gi09 = new TGraph( n09, &vb09[0], &ls09[0] );
  gi09->SetLineColor(4);
  gi09->SetMarkerColor(4);
  gi09->SetMarkerStyle(20);
  gi09->SetMarkerSize(1.4);
  gi09->Draw("Pc"); // without axis option: overlay

  gi21 = new TGraph( n21, &vb21[0], &ls21[0] );
  gi21->SetMarkerColor(414); // green
  gi21->SetMarkerStyle(22);
  gi21->SetMarkerSize(1.6);
  gi21->Draw("Pc"); // without axis option: overlay

  gi42 = new TGraph( n42, &vb42[0], &ls42[0] );
  gi42->SetMarkerColor(419); // dark green
  gi42->SetMarkerStyle(21);
  gi42->SetMarkerSize(1.2);
  gi42->Draw("Pc"); // without axis option: overlay

  gi89 = new TGraph( n89, &vb89[0], &ls89[0] );
  gi89->SetLineColor(2);
  gi89->SetMarkerColor(2);
  gi89->SetMarkerStyle(29);
  gi89->SetMarkerSize(1.9);
  gi89->Draw("Pc"); // without axis option: overlay

  TLegend * ll = new TLegend( 0.56, 0.66, 0.95, 0.90 );
  ll->AddEntry( gi89, "589, 2.1#upoint10^{16} n_{Ka}/cm^{2} ", "p" );
  ll->AddEntry( gi42, "542, 7.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  ll->AddEntry( gi21, "521, 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  ll->AddEntry( gi09, "509, 5.2#upoint10^{15} n_{PS}/cm^{2} ", "p" );
  ll->Draw( "same" );

  c1.Print( "bias.pdf" );

  c1.SetGridx(0);
  c1.SetGridy(0);
  gPad->SetLogy(0);

  //----------------------------------------------------------------------------
  // eff vs bias

  he12 = new
    TH1F( "he12",
	  "RD53A 2.7#upoint10^{15 }n_{PS}/cm^{2};sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  he12->SetStats(kFALSE); // no statistics
  he12->SetMinimum( 80);
  he12->SetMaximum(100);
  he12->Draw();

  ge12 = new TGraph( n12, &vb12[0], &ef12[0] );
  ge12->SetLineColor(2);
  ge12->SetMarkerColor(2);
  ge12->SetMarkerStyle(20);
  ge12->SetMarkerSize(1.4);
  ge12->Draw("Pc"); // without axis option: overlay

  ge1234 = new TGraph( n1234, &vb1234[0], &ef1234[0] );
  ge1234->SetLineColor(kRed+2);
  ge1234->SetMarkerColor(kRed+2);
  ge1234->SetMarkerStyle(21);
  ge1234->SetMarkerSize(1.2);
  ge1234->Draw("Pc"); // without axis option: overlay

  TLegend * lgnd12 = new TLegend( 0.5, 0.18, 0.94, 0.32 );
  lgnd12->AddEntry( ge12,   "512 bias dot, turn 0 ", "pl" );
  lgnd12->AddEntry( ge1234, "512 bias dot, turn 34^{o} ", "pl" );
  lgnd12->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // eff vs bias with turn

  het = new
    TH1F( "het",
	  "RD53A with turn;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  het->SetStats(kFALSE); // no statistics
  het->SetMinimum( 80);
  het->SetMaximum(100);
  het->Draw();

  ge1234->Draw("Pc"); // without axis option: overlay

  ge21 = new TGraph( n21, &vb21[0], &ef21[0] );
  ge21->SetLineColor(414); // green
  ge21->SetMarkerColor(414); // green
  ge21->SetMarkerStyle(22);
  ge21->SetMarkerSize(1.6);
  ge21->Draw("Pc"); // without axis option: overlay

  ge11t = new TGraph( n11t, &vb11t[0], &ef11t[0] );
  ge11t->SetLineColor(1);
  ge11t->SetMarkerColor(1);
  ge11t->SetMarkerStyle(23);
  ge11t->SetMarkerSize(1.6);
  ge11t->Draw("Pc"); // without axis option: overlay

  TLegend * lgndt = new TLegend( 0.3, 0.18, 0.94, 0.40 );
  lgndt->AddEntry( ge21, "521 bias dot, 18^{o}, 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "pl" );
  lgndt->AddEntry( ge1234, "512 bias dot, 34^{o}, 5.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndt->AddEntry( ge11t,  "511 bias dot, 30^{o}, 6.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndt->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // CERN PS

  hePS = new
    TH1F( "hePS",
	  "RD53A CERN PS;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  hePS->SetStats(kFALSE); // no statistics
  hePS->SetMinimum( 80);
  hePS->SetMaximum(100);
  hePS->Draw();

  ge0934->SetMarkerStyle(20);
  ge0934->Draw("Pc"); // without axis option: overlay

  ge1234->Draw("Pc"); // without axis option: overlay

  ge11t->Draw("Pc"); // without axis option: overlay

  TLegend * lgndPS = new TLegend( 0.3, 0.18, 0.94, 0.40 );
  lgndPS->AddEntry( ge0934, "509 default, 34^{o}, 5.2#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndPS->AddEntry( ge1234, "512 bias dot, 34^{o}, 5.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndPS->AddEntry( ge11t,  "511 bias dot, 30^{o}, 6.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndPS->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // all

  hea = new
    TH1F( "hea",
	  "all;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  hea->SetStats(kFALSE); // no statistics
  hea->SetMinimum( 80);
  hea->SetMaximum(100);
  hea->Draw();

  ge0934->Draw("Pc"); // without axis option: overlay

  ge21->Draw("Pc"); // without axis option: overlay

  ge1234->Draw("Pc"); // without axis option: overlay

  ge11t->Draw("Pc"); // without axis option: overlay

  TLegend * lgnda = new TLegend( 0.3, 0.18, 0.94, 0.40 );
  lgnda->AddEntry( ge0934, "509 default, 34^{o}, 5.2#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgnda->AddEntry( ge21, "521 bias dot, 18^{o}, 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "pl" );
  lgnda->AddEntry( ge1234, "512 bias dot, 34^{o}, 5.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgnda->AddEntry( ge11t,  "511 bias dot, 30^{o}, 6.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgnda->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // all legend for externals

  hfx = new
    TH1F( "hfx",
	  "all;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  hfx->SetStats(kFALSE); // no statistics
  hfx->SetMinimum( 80);
  hfx->SetMaximum(100);
  hfx->Draw();

  ge0934->Draw("Pc"); // without axis option: overlay

  ge21->Draw("Pc"); // without axis option: overlay

  ge1234->Draw("Pc"); // without axis option: overlay

  ge11t->Draw("Pc"); // without axis option: overlay

  // legend for externals:

  TLegend * lgndx = new TLegend( 0.28, 0.18, 0.94, 0.38 );
  lgndx->AddEntry( ge0934, "100x25 no bias dot, 5.2#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndx->AddEntry( ge21, "50x50 bias dot, 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "pl" );
  lgndx->AddEntry( ge1234, "100x25 bias dot, 5.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndx->AddEntry( ge11t,  "50x50 bias dot, 6.6#upoint10^{15 }n_{PS}/cm^{2} ", "pl" );
  lgndx->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // charge vs bias, 589

  hq89 = new
    TH1F( "hq89",
	  "2.1#upoint10^{16} n (Ka);sensor bias [V];RD53A Lin mean charge [ToT]",
	  82, 0, 820 ); // axis range
  hq89->SetStats(kFALSE); // no statistics
  hq89->SetMinimum( 0);
  hq89->SetMaximum(13);
  hq89->Draw();

  gq89 = new TGraph( n89, &vb89[0], &q089[0] );
  gq89->SetLineColor(2);
  gq89->SetMarkerColor(2);
  gq89->SetMarkerStyle(29);
  gq89->SetMarkerSize(1.9);
  gq89->Draw("Pc"); // without axis option: overlay

  gq89b = new TGraph( n89b-4, &vb89b[0], &q089b[0] );
  gq89b->SetLineColor(634); // dark red
  gq89b->SetMarkerColor(634);
  gq89b->SetMarkerStyle(22);
  gq89b->SetMarkerSize(1.5);
  gq89b->Draw("Pc"); // without axis option: overlay

  gq89r = new TGraph( n89b-4, &vb89b[0], &q089r[0] );
  gq89r->SetLineColor(1);
  gq89r->SetMarkerColor(1);
  gq89r->SetMarkerStyle(26);
  gq89r->SetMarkerSize(1.5);
  gq89r->Draw("Pc"); // without axis option: overlay

  TLegend * lq89 = new TLegend( 0.50, 0.16, 0.94, 0.32 );
  lq89->AddEntry( gq89, "589, Krum_Curr_Lin 29", "p" );
  lq89->AddEntry( gq89b, "589, Krum_Curr_Lin 20", "p" );
  lq89->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // charge vs bias, defaults:

  hq = new
    TH1F( "hq",
	  "RD53A default q;sensor bias [V];RD53A Lin mean charge [ToT]",
	  82, 0, 820 ); // axis range
  hq->SetStats(kFALSE); // no statistics
  hq->SetMinimum( 0);
  hq->SetMaximum(13);
  hq->Draw();

  gq09 = new TGraph( n09, &vb09[0], &q009[0] );
  gq09->SetLineColor(4);
  gq09->SetMarkerColor(4);
  gq09->SetMarkerStyle(20);
  gq09->SetMarkerSize(1.4);
  gq09->Draw("Pc"); // without axis option: overlay

  gq47 = new TGraph( n47, &vb47[0], &q047[0] );
  gq47->SetLineColor(84); // dark green
  gq47->SetMarkerColor(84); // green
  gq47->SetMarkerStyle(22);
  gq47->SetMarkerSize(1.6);
  gq47->Draw("Pc"); // without axis option: overlay

  gq42 = new TGraph( n42, &vb42[0], &q042[0] );
  gq42->SetLineColor(419); // dark green
  gq42->SetMarkerColor(419); // dark green
  gq42->SetMarkerStyle(23);
  gq42->SetMarkerSize(1.6);
  gq42->Draw("Pc"); // without axis option: overlay

  gq63 = new TGraph( n63, &vb63[0], &q063[0] );
  gq63->SetLineColor(1);
  gq63->SetMarkerColor(1);
  gq63->SetMarkerStyle(20);
  gq63->SetMarkerSize(1.4);
  gq63->Draw("Pc"); // without axis option: overlay

  gq64 = new TGraph( n64, &vb64[0], &q064[0] );
  gq64->SetLineColor(633);
  gq64->SetMarkerColor(633);
  gq64->SetMarkerStyle(21);
  gq64->SetMarkerSize(1.2);
  gq64->Draw("Pc"); // without axis option: overlay

  gq65 = new TGraph( n65, &vb65[0], &q065[0] );
  gq65->SetLineColor(797); // orange
  gq65->SetMarkerColor(797);
  gq65->SetMarkerStyle(23);
  gq65->SetMarkerSize(1.6);
  gq65->Draw("Pc"); // without axis option: overlay

  gq78 = new TGraph( n78, &vb78[0], &q078[0] );
  gq78->SetLineColor(kMagenta+1);
  gq78->SetMarkerColor(kMagenta+1);
  gq78->SetMarkerStyle(22);
  gq78->SetMarkerSize(1.6);
  gq78->Draw("Pc"); // without axis option: overlay

  gq77 = new TGraph( n77, &vb77[0], &q077[0] );
  gq77->SetLineColor(922);
  gq77->SetMarkerColor(922);
  gq77->SetMarkerStyle(29);
  gq77->SetMarkerSize(1.9);
  gq77->Draw("Pc"); // without axis option: overlay

  gq89r->SetLineColor(2);
  gq89r->SetMarkerColor(2);
  gq89r->SetMarkerStyle(29);
  gq89r->SetMarkerSize(1.9);
  gq89r->Draw("Pc"); // without axis option: overlay

  TLegend * lgndq = new TLegend( 0.61, 0.16, 0.95, 0.49 );
  lgndq->AddEntry( gq78, "578: 5.0#upoint10^{15} n (Ka) ", "p" );
  lgndq->AddEntry( gq47, "547: 4.4#upoint10^{15} n (Ka) ", "p" );
  lgndq->AddEntry( gq42, "542: 7.4#upoint10^{15} n (Ka) ", "p" );
  lgndq->AddEntry( gq77, "577: 1.0#upoint10^{16} n (Ka) ", "p" );
  lgndq->AddEntry( gq09, "509: 5.2#upoint10^{15} n (PS) ", "p" );
  lgndq->AddEntry( gq65, "565: 5.0#upoint10^{15} n (Bh) ", "p" );
  lgndq->AddEntry( gq64, "564: 5.0#upoint10^{15} n (Bh) ", "p" );
  lgndq->AddEntry( gq63, "563: 1.0#upoint10^{16} n (Bh) ", "p" );
  lgndq->AddEntry( gq89, "589, 2.1#upoint10^{16} n (Ka) ", "p" );
  lgndq->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // eff vs bias, default designs:

  heff = new
    TH1F( "heff",
	  "RD53A default eff;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  heff->SetStats(kFALSE); // no statistics
  heff->SetMinimum( 90);
  heff->SetMaximum(100);
  heff->Draw();

  TLine l99( 0, 99, 820, 99 );
  l99.SetLineStyle(2);
  l99.Draw("same");

  ge09->SetMarkerColor(4);
  ge09->SetLineColor(4);
  ge09->Draw("Pc"); // without axis option: overlay

  ge47 = new TGraph( n47, &vb47[0], &ef47[0] );
  ge47->SetLineColor(84); // dark green
  ge47->SetMarkerColor(84); // green
  ge47->SetMarkerStyle(22);
  ge47->SetMarkerSize(1.6);
  ge47->Draw("Pc"); // without axis option: overlay

  ge42 = new TGraph( n42, &vb42[0], &ef42[0] );
  ge42->SetLineColor(419); // dark green
  ge42->SetMarkerColor(419); // dark green
  ge42->SetMarkerStyle(23);
  ge42->SetMarkerSize(1.6);
  ge42->Draw("Pc"); // without axis option: overlay

  ge78 = new TGraph( n78, &vb78[0], &ef78[0] );
  ge78->SetLineColor(kMagenta+1);
  ge78->SetMarkerColor(kMagenta+1);
  ge78->SetMarkerStyle(22);
  ge78->SetMarkerSize(1.6);
  ge78->Draw("Pc"); // without axis option: overlay

  ge77 = new TGraph( n77, &vb77[0], &ef77[0] );
  ge77->SetLineColor(922); // gray
  ge77->SetMarkerColor(922);
  ge77->SetMarkerStyle(29);
  ge77->SetMarkerSize(1.9);
  ge77->Draw("Pc"); // without axis option: overlay

  ge64 = new TGraph( n64, &vb64[0], &ef64[0] );
  ge64->SetLineColor(633);
  ge64->SetMarkerColor(633);
  ge64->SetMarkerStyle(21);
  ge64->SetMarkerSize(1.2);
  ge64->Draw("Pc"); // without axis option: overlay

  ge65 = new TGraph( n65, &vb65[0], &ef65[0] );
  ge65->SetLineColor(797); // orange
  ge65->SetMarkerColor(797);
  ge65->SetMarkerStyle(23);
  ge65->SetMarkerSize(1.6);
  ge65->Draw("Pc"); // without axis option: overlay

  ge63 = new TGraph( n63, &vb63[0], &ef63[0] );
  ge63->SetLineColor(1);
  ge63->SetMarkerColor(1);
  ge63->SetMarkerStyle(20);
  ge63->SetMarkerSize(1.4);
  ge63->Draw("Pc"); // without axis option: overlay

  ge89b = new TGraph( n89b-4, &vb89b[0], &ef89b[0] );
  ge89b->SetLineColor(2);
  ge89b->SetMarkerColor(2);
  ge89b->SetMarkerStyle(29);
  ge89b->SetMarkerSize(1.9);
  ge89b->Draw("Pc"); // without axis option: overlay

  TLegend * lgndef = new TLegend( 0.61, 0.16, 0.94, 0.60 );
  lgndef->AddEntry( ge78, "578: 5.0#upoint10^{15} n (Ka) ", "p" );
  lgndef->AddEntry( ge09, "509: 5.2#upoint10^{15} n (PS) ", "p" );
  lgndef->AddEntry( ge47, "547: 4.4#upoint10^{15} n (Ka) ", "p" );
  lgndef->AddEntry( ge42, "542: 7.4#upoint10^{15} n (Ka) ", "p" );
  lgndef->AddEntry( ge77, "577: 1.0#upoint10^{16} n (Ka) ", "p" );
  lgndef->AddEntry( ge65, "565: 5.0#upoint10^{15} n (Bh) ", "p" );
  lgndef->AddEntry( ge64, "564: 5.0#upoint10^{15} n (Bh) ", "p" );
  lgndef->AddEntry( ge63, "563: 1.0#upoint10^{16} n (Bh) ", "p" );
  lgndef->AddEntry( ge89b, "589: 2.1#upoint10^{16} n (Ka) ", "p" );
  lgndef->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // 589 eff vs bias

  he89 = new
    TH1F( "he89",
	  "2.1#upoint10^{16} n (Ka);sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  he89->SetStats(kFALSE); // no statistics
  he89->SetMinimum( 90);
  he89->SetMaximum(100);
  he89->Draw();

  l99.Draw("same");

  ge89b->Draw("Pc"); // without axis option: overlay

  ge89 = new TGraph( n89, &vb89[0], &ef89[0] );
  ge89->SetLineColor(634); // dark red
  ge89->SetMarkerColor(634);
  ge89->SetMarkerStyle(22);
  ge89->SetMarkerSize(1.5);
  ge89->Draw("Pc"); // without axis option: overlay

  TLegend * lgnd89 = new TLegend( 0.60, 0.16, 0.94, 0.36 );
  lgnd89->AddEntry( ge89, "589, thr 1170 e", "p" );
  lgnd89->AddEntry( ge89b, "589, thr 1240 e", "p" );
  lgnd89->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // eff vs bias, Birmingham

  hebh = new
    TH1F( "hebh",
	  "RD53A default eff;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  hebh->SetStats(kFALSE); // no statistics
  hebh->SetMinimum( 90);
  hebh->SetMaximum(100);
  hebh->Draw();

  l99.Draw("same");

  ge09->SetMarkerColor(4);
  ge09->SetLineColor(4);
  ge09->Draw("Pc"); // without axis option: overlay

  ge64->Draw("Pc"); // without axis option: overlay

  ge65->Draw("Pc"); // without axis option: overlay

  ge63->Draw("Pc"); // without axis option: overlay

  TLegend * lgndbh = new TLegend( 0.55, 0.16, 0.94, 0.44 );
  lgndbh->AddEntry( ge09, "25-509: 5.2#upoint10^{15} PS ", "p" );
  lgndbh->AddEntry( ge65, "25-565: 5.0#upoint10^{15} Bh ", "p" );
  lgndbh->AddEntry( ge64, "50-564: 5.0#upoint10^{15} Bh ", "p" );
  lgndbh->AddEntry( ge63, "25-563: 1.0#upoint10^{16} Bh ", "p" );
  lgndbh->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // eff vs bias, Ka

  heka = new
    TH1F( "heka",
	  "RD53A default eff;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  heka->SetStats(kFALSE); // no statistics
  heka->SetMinimum( 90);
  heka->SetMaximum(100);
  heka->Draw();

  l99.Draw("same");
  ge47->Draw("Pc"); // without axis option: overlay
  ge42->Draw("Pc"); // without axis option: overlay
  ge78->Draw("Pc"); // without axis option: overlay
  ge77->Draw("Pc"); // without axis option: overlay

  TLegend * lgndka = new TLegend( 0.55, 0.16, 0.94, 0.44 );
  lgndka->AddEntry( ge47, "50-547: 4.4#upoint10^{15} Ka ", "p" );
  lgndka->AddEntry( ge78, "25-578: 5.0#upoint10^{15} Ka ", "p" );
  lgndka->AddEntry( ge42, "25-542: 7.4#upoint10^{15} Ka ", "p" );
  lgndka->AddEntry( ge77, "25-577: 1.0#upoint10^{16} Ka ", "p" );
  lgndka->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // eff vs q, default designs, vertical

  heq = new
    TH1F( "heq",
	  "RD53A default eff;signal [ToT];RD53A Lin efficiency [%]",
	  27, 0, 13.5 ); // axis range
  heq->SetStats(kFALSE); // no statistics
  heq->SetMinimum( 90);
  heq->SetMaximum(100);
  heq->Draw();

  TLine lq99( 0, 99, 12, 99 );
  lq99.SetLineStyle(2);
  lq99.Draw("same");

  geq09 = new TGraph( n09, &q009[0], &ef09[0] );
  geq09->SetLineColor(4);
  geq09->SetMarkerColor(4);
  geq09->SetMarkerStyle(20);
  geq09->SetMarkerSize(1.4);
  geq09->Draw("Pc"); // without axis option: overlay

  geq47 = new TGraph( n47, &q047[0], &ef47[0] );
  geq47->SetLineColor(84); // dark green
  geq47->SetMarkerColor(84); // green
  geq47->SetMarkerStyle(22);
  geq47->SetMarkerSize(1.6);
  geq47->Draw("Pc"); // without axis option: overlay

  geq42 = new TGraph( n42, &q042[0], &ef42[0] );
  geq42->SetLineColor(419); // dark green
  geq42->SetMarkerColor(419); // dark green
  geq42->SetMarkerStyle(23);
  geq42->SetMarkerSize(1.6);
  geq42->Draw("Pc"); // without axis option: overlay

  geq78 = new TGraph( n78, &q078[0], &ef78[0] );
  geq78->SetLineColor(kMagenta+1);
  geq78->SetMarkerColor(kMagenta+1);
  geq78->SetMarkerStyle(22);
  geq78->SetMarkerSize(1.6);
  geq78->Draw("Pc"); // without axis option: overlay

  geq77 = new TGraph( n77, &q077[0], &ef77[0] );
  geq77->SetLineColor(922); // gray
  geq77->SetMarkerColor(922);
  geq77->SetMarkerStyle(29);
  geq77->SetMarkerSize(1.9);
  geq77->Draw("Pc"); // without axis option: overlay

  geq64 = new TGraph( n64, &q064[0], &ef64[0] );
  geq64->SetLineColor(633);
  geq64->SetMarkerColor(633);
  geq64->SetMarkerStyle(21);
  geq64->SetMarkerSize(1.2);
  geq64->Draw("Pc"); // without axis option: overlay

  geq65 = new TGraph( n65, &q065[0], &ef65[0] );
  geq65->SetLineColor(797); // orangeq
  geq65->SetMarkerColor(797);
  geq65->SetMarkerStyle(23);
  geq65->SetMarkerSize(1.6);
  geq65->Draw("Pc"); // without axis option: overlay

  geq63 = new TGraph( n63, &q063[0], &ef63[0] );
  geq63->SetLineColor(1);
  geq63->SetMarkerColor(1);
  geq63->SetMarkerStyle(20);
  geq63->SetMarkerSize(1.4);
  geq63->Draw("Pc"); // without axis option: overlay

  geq89 = new TGraph( n89b-4, &q089r[0], &ef89b[0] );
  geq89->SetLineColor(2);
  geq89->SetMarkerColor(2);
  geq89->SetMarkerStyle(29);
  geq89->SetMarkerSize(1.9);
  geq89->Draw("Pc"); // without axis option: overlay

  TLegend * leq = new TLegend( 0.58, 0.16, 0.94, 0.64 );
  leq->AddEntry( geq09, "509: 5.2#upoint10^{15} n (PS) ", "p" );
  leq->AddEntry( geq65, "565: 5.0#upoint10^{15} n (Bh) ", "p" );
  leq->AddEntry( geq64, "564: 5.0#upoint10^{15} n (Bh) ", "p" );
  leq->AddEntry( geq42, "542: 7.4#upoint10^{15} n (Ka) ", "p" );
  leq->AddEntry( geq78, "578: 5.0#upoint10^{15} n (Ka) ", "p" );
  leq->AddEntry( geq77, "577: 1.0#upoint10^{16} n (Ka) ", "p" );
  leq->AddEntry( geq63, "563: 1.0#upoint10^{16} n (Bh) ", "p" );
  leq->AddEntry( geq47, "547: 4.4#upoint10^{15} n (Ka) ", "p" );
  leq->AddEntry( geq89, "589: 2.1#upoint10^{16} n (Ka) ", "p" );
  leq->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // anneal 509

  ht = new
    TH1F( "ht",
	  "5.2#upoint10^{15} n (PS);sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  ht->SetStats(kFALSE); // no statistics
  ht->SetMinimum( 90);
  ht->SetMaximum(100);
  ht->Draw();

  l99.Draw("same");

  ge09->SetMarkerColor(4);
  ge09->SetLineColor(4);
  ge09->Draw("Pc"); // without axis option: overlay

  ge09a = new TGraph( n09a, &vb09a[0], &ef09a[0] );
  ge09a->SetLineColor(861); // azure
  ge09a->SetMarkerColor(861);
  ge09a->SetMarkerStyle(29);
  ge09a->SetMarkerSize(1.9);
  ge09a->Draw("Pc"); // without axis option: overlay

  TLegend * lgndn = new TLegend( 0.55, 0.20, 0.94, 0.36 );
  lgndn->AddEntry( ge09, "509 Dec 2018 ", "p" );
  lgndn->AddEntry( ge09a, "509 Jun 2020 ", "p" );
  lgndn->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // anneal 542

  h42 = new
    TH1F( "h42",
	  "7.4#upoint10^{15} n (Ka);sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  h42->SetStats(kFALSE); // no statistics
  h42->SetMinimum( 90);
  h42->SetMaximum(100);
  h42->Draw();

  l99.Draw("same");

  ge42->Draw("Pc"); // without axis option: overlay

  ge42a = new TGraph( n42a, &vb42a[0], &ef42a[0] );
  ge42a->SetLineColor(8);
  ge42a->SetMarkerColor(8);
  ge42a->SetMarkerStyle(29);
  ge42a->SetMarkerSize(1.9);
  ge42a->Draw("Pc"); // without axis option: overlay

  TLegend * lgnd42 = new TLegend( 0.55, 0.20, 0.94, 0.36 );
  lgnd42->AddEntry( ge42, "542 Dec 2019 ", "p" );
  lgnd42->AddEntry( ge42a, "542 Jun 2020 ", "p" );
  lgnd42->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // anneal 521

  h21 = new
    TH1F( "h21",
	  "5.3#upoint10^{15} n (Ka);sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  h21->SetStats(kFALSE); // no statistics
  h21->SetMinimum( 90);
  h21->SetMaximum(100);
  h21->Draw();

  l99.Draw("same");

  ge21->Draw("Pc"); // without axis option: overlay

  ge21a = new TGraph( n21a, &vb21a[0], &ef21a[0] );
  ge21a->SetLineColor(417); // green
  ge21a->SetMarkerColor(417);
  ge21a->SetMarkerStyle(29);
  ge21a->SetMarkerSize(1.9);
  ge21a->Draw("Pc"); // without axis option: overlay

  TLegend * lgnd21 = new TLegend( 0.55, 0.20, 0.94, 0.36 );
  lgnd21->AddEntry( ge21, "521 Mar 2019 ", "p" );
  lgnd21->AddEntry( ge21a, "521 Jun 2020 ", "p" );
  lgnd21->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // q0 521

  hq21 = new
    TH1F( "hq21",
	  "5.3#upoint10^{15} n (Ka);sensor bias [V];RD53A Lin mean charge [ToT]",
	  82, 0, 820 ); // axis range
  hq21->SetStats(kFALSE); // no statistics
  hq21->SetMinimum( 0);
  hq21->SetMaximum(14);
  hq21->Draw();

  gq21 = new TGraph( n21, &vb21[0], &q021[0] );
  gq21->SetLineColor(414);
  gq21->SetMarkerColor(414);
  gq21->SetMarkerStyle(22);
  gq21->SetMarkerSize(1.6);
  gq21->Draw("Pc"); // without axis option: overlay

  gq21a = new TGraph( n21a, &vb21a[0], &q021a[0] );
  gq21a->SetLineColor(417); // green
  gq21a->SetMarkerColor(417);
  gq21a->SetMarkerStyle(29);
  gq21a->SetMarkerSize(1.9);
  gq21a->Draw("Pc"); // without axis option: overlay

  TLegend * lgndq21 = new TLegend( 0.55, 0.20, 0.94, 0.36 );
  lgndq21->AddEntry( gq21, "521 Mar 2019 ", "p" );
  lgndq21->AddEntry( gq21a, "521 Jun 2020 ", "p" );
  lgndq21->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // 589 hot vs bias

  gStyle->SetTitleOffset( 2.0, "y" );
  ht89i = new
    TH1F( "ht89i",
	  "2.1#upoint10^{16} n (Ka);sensor bias [V];hot pixels",
	  82, 0, 820 ); // axis range
  ht89i->SetStats(kFALSE); // no statistics
  ht89i->SetMinimum(  0);
  ht89i->SetMaximum(3e3);
  ht89i->Draw();

  ght89 = new TGraph( n89b-4, &vb89b[0], &msk89[0] );
  ght89->SetLineColor(4);
  ght89->SetMarkerColor(4);
  ght89->SetMarkerStyle(29); // star
  ght89->SetMarkerSize(1.9);
  ght89->Draw("Pc"); // without axis option: overlay

  ght89i = new TGraph( 4, &vb89b[n89b-4], &msk89[n89b-4] );
  ght89i->SetLineColor(2);
  ght89i->SetMarkerColor(2);
  ght89i->SetMarkerStyle(20);
  ght89i->SetMarkerSize(1.4);
  ght89i->Draw("Pc"); // without axis option: overlay

  TLegend * l89h = new TLegend( 0.4, 0.5, 0.8, 0.66 );
  l89h->AddEntry( ght89, "589, thr 1240 e", "p" );
  l89h->AddEntry( ght89i, "589, reduced cooling", "p" );
  l89h->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // 589i eff vs bias

  he89i = new
    TH1F( "he89i",
	  "2.1#upoint10^{16} n (Ka);sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  he89i->SetStats(kFALSE); // no statistics
  he89i->SetMinimum( 90);
  he89i->SetMaximum(100);
  he89i->Draw();

  l99.Draw("same");

  ge89b->SetLineColor(4);
  ge89b->SetMarkerColor(4);
  ge89b->Draw("Pc"); // without axis option: overlay

  ge89i = new TGraph( 4, &vb89b[n89b-4], &ef89b[n89b-4] );
  ge89i->SetLineColor(2);
  ge89i->SetMarkerColor(2);
  ge89i->SetMarkerStyle(20);
  ge89i->SetMarkerSize(1.4);
  ge89i->Draw("Pc"); // without axis option: overlay

  TLegend * l89e = new TLegend( 0.56, 0.20, 0.94, 0.32 );
  l89e->AddEntry( ge89b, "589, thr 1240 e", "p" );
  l89e->AddEntry( ge89i, "589, reduced cooling", "p" );
  l89e->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // done:

  c1.Print( "bias.pdf]" ); // ] closes file
  cout << "evince bias.pdf" << endl;

}
