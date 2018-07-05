function [R_ijkl,gamma]=reimann_curveature(g_ij,vector_of_sym_coordinates)
gamma=sym(zeros([size(g_ij),length(g_ij)]));
for i=1:length(g_ij)
    for j=1:length(g_ij)
        for k=1:length(g_ij)
            a=0;
            b=0;
            c=0;
            for l=1:length(g_ij)
                a=diff(g_ij(k,l),vector_of_sym_coordinates(j))+a;
                b=diff(g_ij(l,j),vector_of_sym_coordinates(k))+b;
                c=diff(g_ij(j,k),vector_of_sym_coordinates(l))+c;
                gamma(i,j,k)=.5*g_ij(i,l)*(a+b-c)+gamma(i,j,k);
            end
        end
    end
end
R_ijkl=sym(zeros(length(g_ij),length(g_ij),length(g_ij),length(g_ij)));
for i=1:length(g_ij)
    for j=1:length(g_ij)
        for k=1:length(g_ij)
            for l=1:length(g_ij)
                a=diff(gamma(i,l,j),vector_of_sym_coordinates(k));
                b=diff(gamma(l,j,k),vector_of_sym_coordinates(i));
                c=0;
                d=0;
                for m=1:length(g_ij)
                    c=gamma(i,k,m)*gamma(m,l,j)+c;
                    d=gamma(i,l,m)*gamma(m,k,j)+d;
                end
                R_ijkl(i,j,k,l)=a-b+c-d;
            end
        end
    end
end
end
