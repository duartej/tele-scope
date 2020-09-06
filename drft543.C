
// Daniel Pitzl, Apr 2019

// edge-on ToT vs depth

// root -l edge543.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

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
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "edge543.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.83, 0.3, 0.97, 0.89 );

  // 50 V:

  TFile * f50 = TFile::Open( "edg53_36095.root" );
  dutbcvsd->SetTitle( "RD53 543 (fresh) edge-on" );
  dutbcvsd->SetMinimum(17);
  dutbcvsd->SetMaximum(23);
  //dutbcvsd->GetXaxis()->SetTitle("depth [mm]");
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(4);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 50 V", "l" );

  // 120 V:
  /*
  TFile * f120 = TFile::Open( "edg53_36094.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(1);
  dutbcvsd->Draw();
  lgnd->AddEntry( dutbcvsd, " 120 V", "l" );
  */
  // 40 V:

  TFile * f40 = TFile::Open( "edg53_36096.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(2);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 40 V", "l" );

  // 30 V:

  TFile * f30 = TFile::Open( "edg53_36097.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(8);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 30 V", "l" );

  // 20 V:

  TFile * f20 = TFile::Open( "edg53_36098.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(93); // orange
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 20 V", "l" );

  // 10 V:

  TFile * f10 = TFile::Open( "edg53_36099.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(6);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 10 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "edge543.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "edge543.pdf]" ); // ] closes file
  cout << "evince edge543.pdf" << endl;

}
