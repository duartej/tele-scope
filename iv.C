
// Daniel Pitzl, Jun 2020

// root -l iv.C

//void iv()
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

  c1.Print( "iv.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  string HASH {"#"};

  //----------------------------------------------------------------------------
  // chip 589: 100x25 default, KIT 2E16

  cout << "try to open dat file for c589";
  ifstream Dstream89( "iv589.dat" );
  if( !Dstream89 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> v89;
  vector <double> i89;

  while( Dstream89.good() && ! Dstream89.eof() ) {

    string rl;
    getline( Dstream89, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double vb;
    double ib;

    strm >> vb;
    strm >> ib;

    v89.push_back(vb-0.02*ib); // 0.02 MOhm drop
    i89.push_back(ib);

  } // while lines

  int n89 = v89.size();
  cout << "    " << n89 << " runs for chip 589i" << endl;

  //----------------------------------------------------------------------------
  // bias scan:

  hiv = new
    TH1F( "hiv",
	  "RD53A 2#upoint10^{16 }n_{Ka}/cm^{2};sensor bias [V];leakage current [#muA]",
	  82, 0, 820 ); // axis range
  hiv->SetStats(kFALSE); // no statistics
  hiv->SetMinimum(0);
  hiv->SetMaximum(999);
  hiv->Draw();

  giv = new TGraph( n89-3, &v89[0], &i89[0] );
  giv->SetLineColor(1);
  giv->SetMarkerColor(1);
  giv->SetMarkerStyle(20);
  giv->SetMarkerSize(1.4);
  giv->Draw("P"); // without axis option: overlay

  TLegend * liv = new TLegend( 0.2, 0.8, 0.7, 0.88 );
  liv->AddEntry( giv, "589: 2#upoint10^{16} n_{Ka} ", "p" );
  liv->Draw( "same" );

  c1.Print( "iv.pdf" );

  //----------------------------------------------------------------------------
  // with warmer points:

  hiw = new
    TH1F( "hiw",
	  "RD53A 2#upoint10^{16 }n_{Ka}/cm^{2};sensor bias [V];leakage current [#muA]",
	  82, 0, 820 ); // axis range
  hiw->SetStats(kFALSE); // no statistics
  hiw->SetMinimum(0);
  hiw->SetMaximum(999);
  hiw->Draw();

  giv->Draw("P"); // without axis option: overlay

  giw = new TGraph( 3, &v89[n89-3], &i89[n89-3] );
  giw->SetLineColor(2);
  giw->SetMarkerColor(2);
  giw->SetMarkerStyle(22);
  giw->SetMarkerSize(1.4);
  giw->Draw("P"); // without axis option: overlay

  TLegend * liw = new TLegend( 0.2, 0.72, 0.7, 0.88 );
  liw->AddEntry( giv, "589: 2#upoint10^{16} n_{Ka} ", "p" );
  liw->AddEntry( giw, "589: less cooling ", "p" );
  liw->Draw( "same" );

  c1.Print( "iv.pdf" );

  //----------------------------------------------------------------------------
  // done:

  c1.Print( "iv.pdf]" ); // ] closes file
  cout << "evince iv.pdf" << endl;

}
