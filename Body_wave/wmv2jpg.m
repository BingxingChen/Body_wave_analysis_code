% save .wmv as .jpg
fileName = 'deg45-3HZfilm0000_clip.avi';
 obj = VideoReader(fileName);
numFrames = obj.NumberOfFrames;% number of frames
 STEP = 1;
 for k = 1:STEP:numFrames% ��ȡ����
      frame = read(obj,k);
 %      imshow(frame);%��ʾ֡
      imwrite(frame,strcat(num2str(k),'.jpg'),'jpg');% ����֡
 end
