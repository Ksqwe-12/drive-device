
char LineOrder[3]={'0','1','2'};
char LO[4]={'l','0','r','0'};
char RoundOrder[3]={'0','1','2'};

void setup() {
  Serial.begin(9600);
  pinMode(LED_BUILTIN, OUTPUT);
  
}

void loop() {
  delay(5);
    int Val_Line = analogRead(A1);
    int Val_Round = analogRead(A0);
    //Serial.println(Val_Round);
    if(Val_Line>700){
      LO[1]=LineOrder[1];
    }else if (Val_Line<300){
     LO[1]=LineOrder[2];
    }else{
      LO[1]=LineOrder[0];
    }
    
    if(Val_Round>700){
      LO[3]=RoundOrder[1];
    }else if (Val_Round<300){
      LO[3]=RoundOrder[2];
    }else{
      LO[3]=RoundOrder[0];
    }
 
    Serial.println(LO);
}
