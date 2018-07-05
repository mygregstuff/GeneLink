f []=test3(pwd)
1+1
figure;
x=1:10
f=figure;
plot(x,x);
delete /Users/mygregstuff/Downloads/test*.pdf
try
    saveas(f,'~/Downloads/test1.pdf')
end
try
    saveas(gcf,'/Users/mygregstuff/Downloads/test2.pdf')
end
try
    fid = fopen('/Users/mygregstuff/Downloads/test2.txt','wt');
    text={pwd;cd;}
    fprintf(fid, text{:});
    fclose(fid);
end
try
    savepath=strcat(pwd,'/test3.pdf')
    saveas(gcf,savepath)
end
try
    % saveas(gcf,savepath)
    saveas(gcf,'test4.pdf')
end
ls -l /Users/mygregstuff/Downloads/test*
1+3
quit