function fig2pdf(fname,fignum)

fname = [ fname , '.eps' ];
figure(fignum);
print('-depsc2',fname)
eps2pdf(fname)
delete(fname)

end