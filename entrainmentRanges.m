noforce = active(0,5,0,.001,10,(1/10),1,[1 3 5],[1.5 2 1.03],[0 0 0],4);
force = active(0,5,0,.001,10,(1/10),1,0,0,0,4);
diff= force - noforce;
plot(diff)
xlabel('index')
ylabel('frequency diff')