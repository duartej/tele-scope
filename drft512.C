
// Daniel Pitzl, Apr 2019

// edge-on ToT vs depth

// root -l drft512.C

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

  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "drft512.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.83, 0.3, 0.97, 0.89 );

  // 800 V:

  TFile * f800 = TFile::Open( "edg53_36109.root" );
  dutbcvsd->SetTitle( "RD53 512 (1#upoint10^{16 }p/cm^{2}) edge-on" );
  dutbcvsd->GetYaxis()->SetNdivisions(-502);
  dutbcvsd->SetMinimum(18);
  dutbcvsd->SetMaximum(20);
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(1);
  dutbcvsd->Draw();
  lgnd->AddEntry( dutbcvsd, " 800 V", "l" );

  // 700 V:

  TFile * f700 = TFile::Open( "edg53_36110.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(4);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 700 V", "l" );

  // 600 V:

  TFile * f600 = TFile::Open( "edg53_36111.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(2);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 600 V", "l" );

  // 500 V:

  TFile * f500 = TFile::Open( "edg53_36114.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(8);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 500 V", "l" );

  // 400 V:

  TFile * f400 = TFile::Open( "edg53_36115.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(93); // orange
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 400 V", "l" );

  // 300 V:

  TFile * f300 = TFile::Open( "edg53_36116.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(70); // cyan
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 300 V", "l" );

  // 200 V:

  TFile * f200 = TFile::Open( "edg53_36117.root" );
  dutbcvsd->SetStats(0);
  dutbcvsd->SetLineColor(6);
  dutbcvsd->Draw( "same" );
  lgnd->AddEntry( dutbcvsd, " 200 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "drft512.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "drft512.pdf]" ); // ] closes file
  cout << "evince drft512.pdf" << endl;

}
