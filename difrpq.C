
// D. Pitzl, Nov 2018
// correlation coefficient profile for scope53

{
  difpqvsym->Clone("difrpqvsym");
  difrpqvsym->Reset();
  difrpqvsym->SetTitle("DIF previous row correlation");
  difrpqvsym->GetYaxis()->SetTitle("DIF previous row signal correlation");
  for( int i = 1; i <= difpqvsym->GetNbinsX(); ++i ) {
    double r = difpqvsym->GetBinContent(i) / sqrt( difppvsym->GetBinContent(i) * difqqvsym->GetBinContent(i) );
    difrpqvsym->Fill( difpqvsym->GetBinCenter(i), r );
    cout << i << " "
	 << difpqvsym->GetBinContent(i) << " "
	 << difppvsym->GetBinContent(i) << " "
	 << difqqvsym->GetBinContent(i) << " "
	 << r << endl;
  }
  difrpqvsym->SetMarkerStyle(20);
  difrpqvsym->Draw("p");
}
