
// Daniel Pitzl, Dec 2018
// RD53A Syyyyyyyn eff vs aut-zero frequency
// root -l auto.C

//------------------------------------------------------------------------------
//void auto()
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleOffset( 1.5, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //                topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );

  c1.Print( "auto.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // RD53A 509 Syn:

  const Int_t n = 5;
  // run            35168 35177 35178 35179 35180
  Float_t frq[n] = {   80,  160,  400,  240,  320 }; // [us] autozero
  Float_t eff[n] = { 92.7, 95.7, 93.2, 96.8, 95.6 }; // 

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // eff vs freq:

  heff = new 
    TH1F( "heff",
	  "efficiency;Sync auto-zero gap [#mus];efficiency [%]",
	  42, 0, 420 ); // axis range
  heff->SetStats(kFALSE); // no statistics
  heff->SetMinimum( 80);
  heff->SetMaximum(100);
  heff->Draw();

  geff = new TGraph( n, frq, eff );
  geff->SetMarkerColor(2);
  geff->SetMarkerStyle(20);
  geff->SetMarkerSize(1.5);
  geff->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.4, 0.2, 0.93, 0.3 );
  lgnd->AddEntry( geff, "RD53A Sync, 8#upoint10^{15} p/cm^{2} ", "p" );
  lgnd->Draw( "same" );

  c1.Print( "auto.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "auto.ps]" ); // ] closes file
  int ierr;
  ierr = system("ps2pdf auto.ps");
  ierr = system("rm -f  auto.ps");
  cout << "acroread auto.pdf" << endl;

}
