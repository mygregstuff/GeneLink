function [R_ijkl,gamma]=reimann_curveature(g_ij,vector_of_sym_coordinates)
gamma=sym(zeros([size(g_ij),length(g_ij)]));
gij=g_ij^-1;
for i=1:length(g_ij)
    for j=1:length(g_ij)
        for k=1:length(g_ij)
            for l=1:length(g_ij)
                a=diff(g_ij(k,l),vector_of_sym_coordinates(j));
                b=diff(g_ij(l,j),vector_of_sym_coordinates(k));
                c=diff(g_ij(j,k),vector_of_sym_coordinates(l));
                gamma(i,j,k)=.5*gij(i,l)*(a+b-c)+gamma(i,j,k);
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
                for m=1:length(g_ij)
                    c=gamma(i,k,m)*gamma(m,l,j);
                    d=gamma(i,l,m)*gamma(m,k,j);
                    R_ijkl(i,j,k,l)=a-b+c-d+R_ijkl(i,j,k,l);
                end
            end
        end
    end
end
end
