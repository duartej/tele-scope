
// Daniel Pitzl, Jun 2020
// bias scans: thr

// root -l scans.C

//void scans()
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

  //--------------------------------------------------------------------
  // square canvas:
  //            topleft x, y, width x, y
  TCanvas c1("c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "scans.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  string HASH {"#"};

  cout << "try to open scans.dat";
  ifstream Dstream( "scans.dat" );
  if( !Dstream ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  TH1I * hthr = new TH1I( "hthr",
			  "RD53A bias scans;threshold [e];bias scans",
			  16, 0, 1600 ); // axis range

  vector <int> vchp;
  vector <int> vthr;

  // Read file by lines:

  while( Dstream.good() && ! Dstream.eof() ) {

    string rl;
    getline( Dstream, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int chip;
    int tag;
    strm >> tag;
    if( tag < 2018 ) {
      chip = tag;
      continue;
    }

    int year = tag;
    int mo;
    strm >> mo;

    int thr;
    strm >> thr;

    int run;
    strm >> run;

    hthr->Fill( thr );
    vthr.push_back(thr);
    vchp.push_back(chip);

  } // while lines

  int nthr = vthr.size();
  cout << "    " << nthr << " scans" << endl;

  //----------------------------------------------------------------------------

  gStyle->SetOptStat(10);
  gStyle->SetStatX(0.6);
  hthr->SetNdivisions(208);
  hthr->Draw();
  c1.Print( "scans.pdf" );

  //----------------------------------------------------------------------------
  // thr vs chip

  TH1I * hchp = new TH1I( "hchp",
			  "RD53A bias scans;single chip module;threshold [e]",
			  99, 500.5, 599.5 ); // axis range
  hchp->SetStats(kFALSE); // no statistics
  hchp->SetMinimum( 700);
  hchp->SetMaximum(1500);
  hchp->GetYaxis()->SetTitleOffset(2.3);
  gPad->SetGridx();
  hchp->Draw();

  TGraph * gchp = new TGraph( nthr, &vchp[0], &vthr[0] );
  gchp->SetLineColor(4);
  gchp->SetMarkerColor(4);
  gchp->SetMarkerStyle(20);
  gchp->SetMarkerSize(1.4);
  gchp->Draw("P"); // without axis option: overlay

  TLegend * lchp = new TLegend( 0.2, 0.2, 0.7, 0.28 );
  lchp->AddEntry( gchp, "RD53A bias scans", "p" );
  //lchp->Draw( "same" );

  c1.Print( "scans.pdf" );

  //----------------------------------------------------------------------------
  // done:

  c1.Print( "scans.pdf]" ); // ] closes file
  cout << "evince scans.pdf" << endl;

}
