 I=imread(strcat(num2str(119),'.jpg'));%��ȡһ֡ͼ��
%  Background=imread('1761.jpg');
%     I=imsubtract(I,Background);
%     I=imabsdiff(I,Background);
%     K=imadd(imabsdiff(I,Background),50);
    imshow(I)
%     return
%     I = rgb2gray(K);%ת��Ϊ�Ҷ�ͼ��
% %     BW = imbinarize(I);
%     BW = ~im2bw(I,Threshhold);%ͼ���ֵ������ 
%     %��ʾͼ��
%     imshow(BW)
    %�����ѡȡĿ������
    h = imrect;
    position = round(wait(h));%[X,Y,W,H]
   
    p1 = position(1);
    p2 = position(2);
    p3 = position(1)+position(3);
    p4 = position(2)+position(4);
    p=[ p1   p2  p3   p4 ]
    
% % p0=p;
distance=sqrt((p0(1)-p(1))^2+(p0(2)-p(2))^2)
velocity=distance/(num_end-num_orign)*30/fish_length0