
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
  gStyle->SetTitleOffset( 2.0, "y" );

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

  cout << "try to open dat file for c509i";
  ifstream Dstream09( "bias509i0.dat" ); // turn 0 = vertical
  if( !Dstream09 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vb09;
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

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb09.push_back(vb);
    ef09.push_back(eff);
    ls09.push_back(100-eff);

  } // while lines

  int n09 = vb09.size();
  cout << "    " << n09 << " runs for chip 509i" << endl;

  //----------------------------------------------------------------------------
  // chip 509 = FTH150P RD53A 100x25 default  5E15 n_eq at CERN PS with turn

  cout << "try to open dat file for c509i34";
  ifstream Dstream0934( "bias509i34.dat" ); // turn 34
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
  // chip 512 = FTH150P RD53A 100x25 bias dot  5.6E15 n_eq at CERN PS

  cout << "try to open dat file for c512i";
  ifstream Dstream12( "bias512i.dat" ); // turn 0 = vertical
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
  // chip 521 = FTH150P RD53A 100x25 bias dot wiggle  5E15 n_eq at KaZ

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

    strm >> vb;
    strm >> run;
    strm >> eff;

    vb21.push_back(vb);
    ef21.push_back(eff);
    ls21.push_back(100-eff); // loss = ineff

  } // while lines

  int n21 = vb21.size();
  cout << "    " << n21 << " runs for chip 521i" << endl;

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
  ge09->SetMarkerColor(2);
  ge09->SetLineColor(2);
  ge09->SetMarkerStyle(20);
  ge09->SetMarkerSize(1.5);
  ge09->Draw("Pc"); // without axis option: overlay

  ge0934 = new TGraph( n0934, &vb0934[0], &ef0934[0] );
  ge0934->SetLineColor(4);
  ge0934->SetMarkerColor(4);
  ge0934->SetMarkerStyle(24);
  ge0934->SetMarkerSize(1.5);
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
  hi->SetTitleOffset( 1.7, "y" );
  hi->Draw();

  gi09 = new TGraph( n09, &vb09[0], &ls09[0] );
  gi09->SetMarkerColor(4);
  gi09->SetMarkerStyle(20);
  gi09->SetMarkerSize(1.5);
  gi09->Draw("Pc"); // without axis option: overlay

  gi21 = new TGraph( n21, &vb21[0], &ls21[0] );
  gi21->SetMarkerColor(413); // green
  gi21->SetMarkerStyle(22);
  gi21->SetMarkerSize(1.5);
  gi21->Draw("Pc"); // without axis option: overlay

  gi42 = new TGraph( n42, &vb42[0], &ls42[0] );
  gi42->SetMarkerColor(419); // dark green
  gi42->SetMarkerStyle(21);
  gi42->SetMarkerSize(1.5);
  gi42->Draw("Pc"); // without axis option: overlay

  TLegend * ll = new TLegend( 0.34, 0.77, 0.94, 0.89 );
  ll->AddEntry( gi09, "509, 5.2#upoint10^{15} n_{PS}/cm^{2} ", "p" );
  ll->AddEntry( gi21, "521, 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  ll->AddEntry( gi42, "542, 7.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
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
  ge12->SetMarkerSize(1.5);
  ge12->Draw("Pc"); // without axis option: overlay

  ge1234 = new TGraph( n1234, &vb1234[0], &ef1234[0] );
  ge1234->SetLineColor(kRed+2);
  ge1234->SetMarkerColor(kRed+2);
  ge1234->SetMarkerStyle(21);
  ge1234->SetMarkerSize(1.5);
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
  ge21->SetLineColor(413); // green
  ge21->SetMarkerColor(413); // green
  ge21->SetMarkerStyle(22);
  ge21->SetMarkerSize(1.5);
  ge21->Draw("Pc"); // without axis option: overlay

  ge11t = new TGraph( n11t, &vb11t[0], &ef11t[0] );
  ge11t->SetLineColor(1);
  ge11t->SetMarkerColor(1);
  ge11t->SetMarkerStyle(23);
  ge11t->SetMarkerSize(1.5);
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
  // Ka

  heKa = new
    TH1F( "heKa",
	  "RD53A KaZ;sensor bias [V];RD53A Lin efficiency [%]",
	  82, 0, 820 ); // axis range
  heKa->SetStats(kFALSE); // no statistics
  heKa->SetMinimum( 90);
  heKa->SetMaximum(100);
  heKa->Draw();

  ge21->Draw("Pc"); // without axis option: overlay

  ge47 = new TGraph( n47, &vb47[0], &ef47[0] );
  ge47->SetLineColor(84); // dark green
  ge47->SetMarkerColor(84); // green
  ge47->SetMarkerStyle(22);
  ge47->SetMarkerSize(1.5);
  ge47->Draw("Pc"); // without axis option: overlay

  ge42 = new TGraph( n42, &vb42[0], &ef42[0] );
  ge42->SetLineColor(419); // dark green
  ge42->SetMarkerColor(419); // dark green
  ge42->SetMarkerStyle(23);
  ge42->SetMarkerSize(1.5);
  ge42->Draw("Pc"); // without axis option: overlay

  TLegend * lgndKa = new TLegend( 0.4, 0.18, 0.94, 0.40 );
  lgndKa->AddEntry( ge47, "547: 4.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  lgndKa->AddEntry( ge21, "521: 5.3#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  lgndKa->AddEntry( ge42, "542: 7.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  lgndKa->Draw( "same" );

  c1.Print( "bias.pdf" );

  //----------------------------------------------------------------------------
  // Ka charge

  hqKa = new
    TH1F( "hqKa",
	  "RD53A KaZ;sensor bias [V];RD53A Lin mean charge [ToT]",
	  82, 0, 820 ); // axis range
  hqKa->SetStats(kFALSE); // no statistics
  hqKa->SetMinimum( 0);
  hqKa->SetMaximum(15);
  hqKa->Draw();

  gq47 = new TGraph( n47, &vb47[0], &q047[0] );
  gq47->SetLineColor(84); // dark green
  gq47->SetMarkerColor(84); // green
  gq47->SetMarkerStyle(22);
  gq47->SetMarkerSize(1.5);
  gq47->Draw("Pc"); // without axis option: overlay

  gq42 = new TGraph( n42, &vb42[0], &q042[0] );
  gq42->SetLineColor(419); // dark green
  gq42->SetMarkerColor(419); // dark green
  gq42->SetMarkerStyle(23);
  gq42->SetMarkerSize(1.5);
  gq42->Draw("Pc"); // without axis option: overlay

  TLegend * lgdqKa = new TLegend( 0.4, 0.18, 0.94, 0.34 );
  lgdqKa->AddEntry( ge47, "547: 4.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  lgdqKa->AddEntry( ge42, "542: 7.4#upoint10^{15} n_{Ka}/cm^{2} ", "p" );
  lgdqKa->Draw( "same" );

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
  // done:

  c1.Print( "bias.pdf]" ); // ] closes file
  cout << "evince bias.pdf" << endl;

}
