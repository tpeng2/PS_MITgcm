for i=1:length(texth)
textstr=texth(i);
textstrnew=num2str(textstr,'%1.1g');
set(texth(i),'String',textstrnew);
end
%% 
levellist=h.LevelList;
levellist=levellist(levellist>=5e-4);
h.LevelList=levellist;