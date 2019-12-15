
// D. Pitzl, Nov 2018
// correlation coefficient profile for scope53

{
  linpqvsym->Clone("linrpqvsym");
  linrpqvsym->Reset();
  linrpqvsym->SetTitle("LIN previous row correlation");
  linrpqvsym->GetYaxis()->SetTitle("LIN previous row signal correlation");
  for( int i = 1; i <= linpqvsym->GetNbinsX(); ++i ) {
    double r = linpqvsym->GetBinContent(i) / sqrt( linppvsym->GetBinContent(i) * linqqvsym->GetBinContent(i) );
    linrpqvsym->Fill( linpqvsym->GetBinCenter(i), r );
    cout << i << " "
	 << linpqvsym->GetBinContent(i) << " "
	 << linppvsym->GetBinContent(i) << " "
	 << linqqvsym->GetBinContent(i) << " "
	 << r << endl;
  }
  linrpqvsym->SetMarkerStyle(20);
  linrpqvsym->Draw("p");
}
