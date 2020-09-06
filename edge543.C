
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
  dutqvsd->SetTitle( "RD53 543 (fresh) edge-on" );
  dutqvsd->SetMaximum(13);
  dutqvsd->GetXaxis()->SetTitle("depth [mm]");
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(4);
  dutqvsd->Draw( "same" );
  lgnd->AddEntry( dutqvsd, " 50 V", "l" );

  // 120 V:
  /*
  TFile * f120 = TFile::Open( "edg53_36094.root" );
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(1);
  dutqvsd->Draw();
  lgnd->AddEntry( dutqvsd, " 120 V", "l" );
  */
  // 40 V:

  TFile * f40 = TFile::Open( "edg53_36096.root" );
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(2);
  dutqvsd->Draw( "same" );
  lgnd->AddEntry( dutqvsd, " 40 V", "l" );

  // 30 V:

  TFile * f30 = TFile::Open( "edg53_36097.root" );
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(8);
  dutqvsd->Draw( "same" );
  lgnd->AddEntry( dutqvsd, " 30 V", "l" );

  // 20 V:

  TFile * f20 = TFile::Open( "edg53_36098.root" );
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(93); // orange
  dutqvsd->Draw( "same" );
  lgnd->AddEntry( dutqvsd, " 20 V", "l" );

  // 10 V:

  TFile * f10 = TFile::Open( "edg53_36099.root" );
  dutqvsd->SetStats(0);
  dutqvsd->SetLineColor(6);
  dutqvsd->Draw( "same" );
  lgnd->AddEntry( dutqvsd, " 10 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "edge543.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "edge543.pdf]" ); // ] closes file
  cout << "evince edge543.pdf" << endl;

}
