//#include<SoftwareSerial.h>
//夹爪使能
// #define enVa 12
// #define enVb 13

class FuzzyPID
{
public:
    FuzzyPID();
    ~FuzzyPID();
    void Get_grad_membership(float erro, float erro_c);
    float Quantization(float maximum, float minimum, float x);
    float Inverse_quantization(float maximum, float minimum, float qvalues);
    void GetSumGrad();
    void GetOUT();

    int changeOut(double out,double sss);

    float FuzzyPIDcontroller(float e_max, float e_min, float ec_max, float ec_min, float kp_max, float kp_min, float erro, float erro_c, float ki_max, float ki_min,float kd_max, float kd_min,float erro_pre, float errp_ppre);
    const int  num_area = 8; //划分区域个数

    float e_membership_values[7] = {-3,-2,-1,0,1,2,3}; //输入e的隶属值
    float ec_membership_values[7] = { -3,-2,-1,0,1,2,3 };//输入de/dt的隶属值
    float kp_menbership_values[7] = { -3,-2,-1,0,1,2,3 };//输出增量kp的隶属值
    float ki_menbership_values[7] = { -3,-2,-1,0,1,2,3 }; //输出增量ki的隶属值
    float kd_menbership_values[7] = { -3,-2,-1,0,1,2,3 };  //输出增量kd的隶属值
    float fuzzyoutput_menbership_values[7] = { -3,-2,-1,0,1,2,3 };

    //int menbership_values[7] = {-3,-};
    float kp;                       //PID参数kp
    float ki;                       //PID参数ki
    float kd;                       //PID参数kd
    float qdetail_kp;               //增量kp对应论域中的值
    float qdetail_ki;               //增量ki对应论域中的值
    float qdetail_kd;               //增量kd对应论域中的值
    float qfuzzy_output;
    float detail_kp;                //输出增量kp
    float detail_ki;                //输出增量ki
    float detail_kd;                //输出增量kd
    float fuzzy_output;
    float qerro;                    //输入e对应论域中的值
    float qerro_c;                  //输入de/dt对应论域中的值
    float errosum;
    float e_gradmembership[2];      //输入e的隶属度
    float ec_gradmembership[2];     //输入de/dt的隶属度
    int e_grad_index[2];            //输入e隶属度在规则表的索引
    int ec_grad_index[2];           //输入de/dt隶属度在规则表的索引
    float gradSums[7] = {0,0,0,0,0,0,0};
    float KpgradSums[7] = { 0,0,0,0,0,0,0 };   //输出增量kp总的隶属度
    float KigradSums[7] = { 0,0,0,0,0,0,0 };   //输出增量ki总的隶属度
    float KdgradSums[7] = { 0,0,0,0,0,0,0 };   //输出增量kd总的隶属度
    int NB = -3, NM = -2, NS = -1, ZO = 0, PS = 1, PM = 2, PB = 3; //论域隶属值

    int  Kp_rule_list[7][7] = { {PB,PB,PM,PM,PS,ZO,ZO},     //kp规则表
                                {PB,PB,PM,PS,PS,ZO,NS},
                                {PM,PM,PM,PS,ZO,NS,NS},
                                {PM,PM,PS,ZO,NS,NM,NM},
                                {PS,PS,ZO,NS,NS,NM,NM},
                                {PS,ZO,NS,NM,NM,NM,NB},
                                {ZO,ZO,NM,NM,NM,NB,NB} };

    int  Ki_rule_list[7][7] = { {NB,NB,NM,NM,NS,ZO,ZO},     //ki规则表
                                {NB,NB,NM,NS,NS,ZO,ZO},
                                {NB,NM,NS,NS,ZO,PS,PS},
                                {NM,NM,NS,ZO,PS,PM,PM},
                                {NM,NS,ZO,PS,PS,PM,PB},
                                {ZO,ZO,PS,PS,PM,PB,PB},
                                {ZO,ZO,PS,PM,PM,PB,PB} };

    int  Kd_rule_list[7][7] = { {PS,NS,NB,NB,NB,NM,PS},    //kd规则表
                                {PS,NS,NB,NM,NM,NS,ZO},
                                {ZO,NS,NM,NM,NS,NS,ZO},
                                {ZO,NS,NS,NS,NS,NS,ZO},
                                {ZO,ZO,ZO,ZO,ZO,ZO,ZO},
                                {PB,NS,PS,PS,PS,PS,PB},
                                {PB,PM,PM,PM,PS,PS,PB} };


};

FuzzyPID::FuzzyPID()  //构造函数
{
    kp = 0;
    ki = 0;
    kd = 0;
    fuzzy_output = 0;
    qdetail_kp = 0;
    qdetail_ki = 0;
    qdetail_kd = 0;
    qfuzzy_output = 0;
    errosum = 0;
}

FuzzyPID::~FuzzyPID()//析构函数
{
}

//输入e与de/dt隶属度计算函数///
void FuzzyPID::Get_grad_membership(float erro,float erro_c)
{
    if (erro > e_membership_values[0] && erro < e_membership_values[6])
    {
        for (int i = 0; i < num_area - 2; i++)
        {
            if (erro >= e_membership_values[i] && erro <= e_membership_values[i + 1])
            {

                //隶属度计算
                e_gradmembership[0] = -(erro - e_membership_values[i + 1]) / (e_membership_values[i + 1] - e_membership_values[i]);
                e_gradmembership[1] = 1+(erro - e_membership_values[i + 1]) / (e_membership_values[i + 1] - e_membership_values[i]);

                //所属区间
                e_grad_index[0] = i;
                e_grad_index[1] = i + 1;
                break;
            }
        }
    }
    else
    {
        if (erro <= e_membership_values[0])
        {
            e_gradmembership[0] = 1;
            e_gradmembership[1] = 0;
            e_grad_index[0] = 0;
            e_grad_index[1] = -1;
        }
        else if (erro >= e_membership_values[6])
        {
            e_gradmembership[0] = 1;
            e_gradmembership[1] = 0;
            e_grad_index[0] = 6;
            e_grad_index[1] = -1;
        }
    }

    if (erro_c > ec_membership_values[0] && erro_c < ec_membership_values[6])
    {
        for (int i = 0; i < num_area - 2; i++)
        {
            if (erro_c >= ec_membership_values[i] && erro_c <= ec_membership_values[i + 1])
            {
                ec_gradmembership[0] = -(erro_c - ec_membership_values[i + 1]) / (ec_membership_values[i + 1] - ec_membership_values[i]);
                ec_gradmembership[1] = 1 + (erro_c - ec_membership_values[i + 1]) / (ec_membership_values[i + 1] - ec_membership_values[i]);
                ec_grad_index[0] = i;
                ec_grad_index[1] = i + 1;
                break;
            }
        }
    }
    else
    {
        if (erro_c <= ec_membership_values[0])
        {
            ec_gradmembership[0] = 1;
            ec_gradmembership[1] = 0;
            ec_grad_index[0] = 0;
            ec_grad_index[1] = -1;
        }
        else if (erro_c >= ec_membership_values[6])
        {
            ec_gradmembership[0] = 1;
            ec_gradmembership[1] = 0;
            ec_grad_index[0] = 6;
            ec_grad_index[1] = -1;
        }
    }

}

//获取输出增量kp,ki,kd的总隶属度/
void FuzzyPID::GetSumGrad()
{
    //连续计算时清空
    for (int i = 0; i <7; i++)
    {
        KpgradSums[i] = 0;
        KigradSums[i] = 0;
        KdgradSums[i] = 0;

    }

   //两两组合，四种情况相累计
  for (int i=0;i<2;i++)
  {
      if (e_grad_index[i] == -1)
      {
       continue;
      }
      for (int j = 0; j < 2; j++)
      {
          if (ec_grad_index[j] != -1)
          {
              //规则从属
              int indexKp = Kp_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;
              int indexKi = Ki_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;
              int indexKd = Kd_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;
              //gradSums[index] = gradSums[index] + (e_gradmembership[i] * ec_gradmembership[j])* Kp_rule_list[e_grad_index[i]][ec_grad_index[j]];

              //根据从属计算隶属度
              KpgradSums[indexKp]= KpgradSums[indexKp] + (e_gradmembership[i] * ec_gradmembership[j]);
              KigradSums[indexKi] = KigradSums[indexKi] + (e_gradmembership[i] * ec_gradmembership[j]);
              KdgradSums[indexKd] = KdgradSums[indexKd] + (e_gradmembership[i] * ec_gradmembership[j]);
          }
          else
          {
            continue;
          }

      }
  }

}

//计算输出增量kp,kd,ki对应论域值//
void FuzzyPID::GetOUT()
{
    for (int i = 0; i < num_area - 1; i++)
    {
        qdetail_kp += kp_menbership_values[i] * KpgradSums[i];
        qdetail_ki += ki_menbership_values[i] * KigradSums[i];
        qdetail_kd+= kd_menbership_values[i] * KdgradSums[i];
    }
}

int FuzzyPID::changeOut(double out,double sss)
{
    int times=0;
    int a=0;
    if(sss<0)
    {
        a=1;
    }else {
        a=0;
    }
    if(out>0)
    {
        times = out/0.53+a;
    }else{
        times = out/0.1-1;
    }
    return times;
}


//模糊PID控制实现函数/
float FuzzyPID::FuzzyPIDcontroller(float e_max, float e_min, float ec_max, float ec_min, float kp_max, float kp_min, float erro, float erro_c,float ki_max,float ki_min,float kd_max,float kd_min,float erro_pre,float errp_ppre)
{
    errosum += erro;
    qerro = Quantization(e_max, e_min, erro);
    qerro_c = Quantization(ec_max, ec_min, erro_c);
    Get_grad_membership(qerro, qerro_c);
    GetSumGrad();
    GetOUT();
    detail_kp = Inverse_quantization(kp_max, kp_min, qdetail_kp);
    detail_ki = Inverse_quantization(ki_max, ki_min, qdetail_ki);
    detail_kd = Inverse_quantization(kd_max, kd_min, qdetail_kd);
    qdetail_kd = 0;
    qdetail_ki = 0;
    qdetail_kp = 0;
    kp = kp + detail_kp;
    ki = ki + detail_ki;
    kd = kd + detail_kd;
  detail_kp = 0;
  detail_ki=0;
  detail_kd=0;
  //qDebug()<<"error:"<<erro<<"\t"<<"erro_pre:"<<erro_pre;
  double sss = erro;
  erro=abs(erro);
  erro_pre=abs(erro_pre);
  errp_ppre= abs(errp_ppre);
  //qDebug()<<"erro:"<<erro<<"\t"<<"erro_pre:"<<erro_pre<<"\t"<<"errp_ppre:"<<errp_ppre;
  float output = kp * erro + ki * (erro - 2 * erro_pre + errp_ppre) + kd*(erro - erro_pre);
  output=changeOut(output,sss);
  //qDebug()<<"kp:"<<kp<<"\t"<<"ki:"<<ki<<"\t"<<"kd:"<<kd<<"out:"<<output;
  return output;
}

///区间映射函数///
float FuzzyPID::Quantization(float maximum,float minimum,float x)
{
    float qvalues= 6.0 *(x-minimum)/(maximum - minimum)-3;
    return qvalues;
}

//反区间映射函数
float FuzzyPID::Inverse_quantization(float maximum, float minimum, float qvalues)
{
    float x = (maximum - minimum) *(qvalues + 3)/6 + minimum;
    return x;
}
































// //夹爪完成标志
// #define isCalA A0
// #define isCalB A1

// //夹爪命令
// #define gA A2
// #define gB A3

//电机使能
 #define enPin 7

 //夹爪命令
 #define A_chong A0
 #define A_pai A1
 #define B_chong A2
 #define B_pai A3

int tiao_flag = 0;

//丝杆电机每次步数
int STEPS_PER_REV = 200;
const int STEPS_PER_REV_1 =200;//长距离
const int STEPS_PER_REV_2 =100;//短距离
//丝杆速度
int X_speed = 750;
const int X_speed_1 = 750;//快速
const int X_speed_2 = 1500;//慢速
//反馈步进电机速度
const int Y_speed = 2000;
//丝杆引脚
const int dirPin = 2;
const int stepPin = 3;
//反馈步进电机引脚
const int dirPin_2 = 4;
const int stepPin_2 = 5;
//直线指令
const int order_go = 1;//前进
const int order_back = -1;//后退
const int order_stop_line = 0;//静止
//旋转指令
const int order_clock = 1;//顺时针旋转
const int order_counter = -1;//逆时针旋转
const int order_stop_round = 0;//静止
//用于接收指令
int Line_Order = 0;
int Round_Order = 0;
//标志位
int gripA=0;//0为没有夹持，1为夹持
int gripB=0;
int wanbi =0;//命令完整接收
//int grip = 0;//夹取标志位
int line_last=0;//上一时刻的直线命令
int round_last=0;//上一时刻的旋转命令

//相对位置
const int max_dis = 7000;//最大限程
const int max_rad = 2500;//最大旋转角度
int line_current_pos = 0;//当前直线位置
int round_current_pos = 0;//当前周线位置

//命令包
float myorder[7];//CurrentActual1 CurrentActual2 CurrentTarget1 CurrentTarget2 LineOrder RoundOrder isAuto
int isAuto = 0;//是否自动控制


 FuzzyPID *myfuzzypid;
 FuzzyPID *myfuzzypid2;

//模糊PID参数
    float Target = 1.5;
    float actual =0;
    float e_max =4;
    float e_min = -4;
    float ec_max = 4;
    float ec_min = -4;

    float kp_max =0.8;
    float kp_min = -0.8;
    float ki_max = 0.1;
    float ki_min = -0.1;
    float kd_max = 0.01;
    float kd_min = -0.01;

    float erro;
    float erro_c;
    float erro_pre = 0;//上一次误差
    float erro_ppre = 0;

int sqq=1;

void setup() {
  //电机引脚初始化
  pinMode(enPin,OUTPUT);
  digitalWrite(enPin, LOW);
  //丝杆引脚
  pinMode(dirPin,OUTPUT);
  pinMode(stepPin,OUTPUT);
  digitalWrite(dirPin, LOW);
  digitalWrite(stepPin, LOW);
  //步进电机引脚
  pinMode(dirPin_2,OUTPUT);
  pinMode(stepPin_2,OUTPUT);
  digitalWrite(dirPin_2, LOW);
  digitalWrite(stepPin_2, LOW);

//   //夹爪使能
//   pinMode(enVa,OUTPUT);
//   digitalWrite(enVa, LOW);//给a使能
//   pinMode(enVb,OUTPUT);
//   digitalWrite(enVb,LOW);
//   //夹爪命令
//   pinMode(gA,OUTPUT);
//   digitalWrite(gA, LOW);//给a充气
//   pinMode(gB,OUTPUT);
//   digitalWrite(gB,LOW);

  //夹爪使能
  pinMode(A_chong,OUTPUT);
  digitalWrite(A_chong, LOW);//给a使能
  pinMode(A_pai,OUTPUT);
  digitalWrite(A_pai,LOW);
  //夹爪命令
  pinMode(B_chong,OUTPUT);
  digitalWrite(B_chong, LOW);//给a充气
  pinMode(B_pai,OUTPUT);
  digitalWrite(B_pai,LOW);
  



  myfuzzypid = new FuzzyPID;//控制夹爪1的模糊pid
  myfuzzypid2 = new FuzzyPID;//控制夹爪2的模糊pid





  //初始化蓝牙 
  Serial.begin(115200);
  pinMode(LED_BUILTIN,OUTPUT);
  //初始化完毕
  digitalWrite(LED_BUILTIN, HIGH);   // turn the LED on (HIGH is the voltage level)
  //delay(5000);
  delay(5000); 
  digitalWrite(LED_BUILTIN, LOW);    // turn the LED off by making the voltage LOW
  
}
 
// the loop function runs over and over again forever
void loop() {
  while(Serial.read() >= 0){}//用于清空缓冲区数据
  delay(50);
  if(Serial.available()>0)
  {
      
      jieshou();//读取命令包
      if (wanbi==1)
      {
          if(sqq==1)
          {
            tiaozheng();
            gripA=1;
            //Serial.println("?");
            sqq=0;
          }
          job();//对所有命令进行响应
          wanbi=0;
      } 
      //Serial.println("?");
  }
}

void jieshou()
{
    while (Serial.available()>0)
    {
      char head = Serial.read();
      if(head=='@')
      {
        fenxi();//读取设定的命令个数
        char tail = Serial.read();
        if(tail=='#'){
              wanbi=1;
              while(Serial.read() >= 0){}
              break;     
              //Serial.println("?");
      } 
      }
    }

}

void fenxi(){
    for(int i=0;i<7;i++)
        {
          myorder[i]=Serial.parseFloat();  
        }
        printpares();//将上数据各自归位分配 
    
}

void job()
{   
      //解析命令
      read_package();//解析直线和周向运动命令
       



    if(!tiao_flag)//判断是否需要调整夹持力
    {   
        //不需要
        Move_Line(Line_Order);//直线运动
        Move_Round(Round_Order);//周向运动
        //状态更新
        line_last = Line_Order;
        round_last = Round_Order;
    }else{
        //需要调整
        if (Target==0)
        {
            if (actual>0.5)//
            {
                open();
                gripA = 0;
            }
        }
        else
        {
            open();//放松之后再调整
           tiaozheng();//使用模糊PID调整夹持力
        }
    }
      
      //状态更新
      line_last = Line_Order;
      round_last = Round_Order;
}

void printpares(){
  //float myorder[6] aa ba at bt lo ro 
   
   tiao_flag = myorder[1];
   Line_Order = myorder[4];
   Round_Order = myorder[5];
   actual = myorder[0];
   Target = myorder[2];
   isAuto = myorder[6];
   if (isAuto==1)//自动与手动的速度不一样，此处进行调整
   {
       STEPS_PER_REV = STEPS_PER_REV_1;
       X_speed = X_speed_1;
   }else
   {
       STEPS_PER_REV = STEPS_PER_REV_2;
       X_speed = X_speed_2;
   } 

}

void tiaozheng()
{
    erro =actual-Target ;
    erro_c = erro - erro_pre;
    float u;
    while(abs(erro)>0.1)//力到达设定范围内即停止调整
    {
        //计算出需要冲/放的次数
        u = myfuzzypid->FuzzyPIDcontroller(e_max, e_min, ec_max, ec_min, kp_max, kp_min, erro, erro_c,ki_max,ki_min,kd_max,kd_min,erro_pre,erro_ppre);
        //Serial.println(u);
        if (u>0)
        {
            while (u--)
            {
                Sachong();
                if(u>5)
                {
                  break;
                  }
            }           
        }else{
            while (abs(u))
            {
                Sapai();
                u++;
            }
        }
        erro_ppre = erro_pre;
        erro_pre = erro;
        while(Serial.read() >= 0){}//清空缓冲区数据
        delay(50);
        jieshou();//再次接收命令包
        //Serial.println("jieshou");
        wanbi=0;
        erro = actual-Target;//本次误差
        erro_c= erro - erro_pre;
    }


}

int Move_Line(int order)
{
    if(order==line_last)//动作指令前后一致
    {
         if(order_stop_line == order)//停止指令
        {
            return 1;
        }
        if(order_go == order)//前进指令
        {
            if(gripA==1)//夹紧
            {
                if((max_dis-line_current_pos)>0)//前进
                {
                    Go_ahead();
                    return 1;
                }
                if((max_dis-line_current_pos)<=0)//到底
                {
                    open();//放松
                    SqueezeB();
                    
                    Go_back();
                    gripA=0;
                    gripB=1;
                    return 1;
                }
            }
            if(gripA==0)
            {
                if(line_current_pos>0)
                {
                    Go_back();
                    return 1;
                }
                if(line_current_pos<=0)
                {
                    openB();
                    Squeeze();//夹紧
                    
                    Go_ahead();
                    gripA=1;
                    gripB=0;
                    return 1;
                }             
            }
        }

        if(order_back==order)//后退指令
        {
             if(gripA==1)//加紧
            {
                if(line_current_pos>0)//后退
                {
                    Go_back();
                    return 1;
                }
                if(line_current_pos<=0)//到底
                {
                  open();//松弛
                    SqueezeB();
                    
                    Go_ahead();
                    gripA=0;
                    gripB=1;
                    return 1;
                }
            }

            if(gripA==0)
            {
                if((max_dis-line_current_pos)>0)
                {
                    Go_ahead();
                    return 1;
                }
                if((max_dis-line_current_pos)<=0)
                {
                  openB();
                    Squeeze();//夹紧
                    
                    Go_back();
                    gripA=1;
                    gripB=0;
                    return 1;
                }             
            }
        }
        
    }


    if(order!=line_last)//动作指令前后不一致
    {
        if(order_stop_line == order)//停止指令
        {
            return 1;
        }

        if(order_go == order)//前进指令
        {
           if(gripA==1)
           {
               if((max_dis-line_current_pos)>0)
               {
                   Go_ahead();
                   return 1;
               }
               if((max_dis-line_current_pos)<=0)
               {
                open();
                   SqueezeB();
                   
                   Go_back();
                   gripA=0;
                   gripB=1;
                   return 1;
               }
           }
           if(gripA==0)
           {
               if((max_dis-line_current_pos)>0)
               {
                openB();
                   Squeeze();
                   
                   Go_ahead();
                   gripA=1;
                   gripB=0;
                   //Go_ahead();
                   return 1;
               }
               if((max_dis-line_current_pos)<=0)
               {
                   Go_back();
                   return 1;
               }
           }
        }

        if(order_back == order)//后退指令
        {
            if(gripA==1)
           {
               if(line_current_pos>0)
               {
                   Go_back();
                   return 1;
               }
               if(line_current_pos<=0)
               {
                open();
                   SqueezeB();
                   
                   Go_ahead();
                   gripA=0;
                   gripB=1;
                   return 1;
               }
           }
           if(gripA==0)
           {
               if(line_current_pos>0)
               {
                openB();
                   Squeeze();
                   
                   //Go_ahead();
                   gripA=1;
                   gripB=0;
                   return 1;
               }
               if(line_current_pos<=0)
               {
                   Go_ahead();
               }
           }
        }
    }
}

int Move_Round(int order)
{

    if(order==round_last)//动作指令前后一致
    {    if(order_stop_round == order)//停止指令
        {
            return 1;
        }

        if(order_clock == order)//前进指令
        {
           if(gripB==1)
           {
               if((max_rad-round_current_pos)>0)
               {
                   Go_Counter();
                   return 1;
               }
               if((max_rad-round_current_pos)<=0)
               { 
                openB();
                   Squeeze();
                   
                   Go_Clock();
                   gripB=0;
                   gripA=1;
                   return 1;
               }
           }
           if(gripB==0)
           {
               if(round_current_pos>0)
               {
                   Go_Clock();
                   return 1;
               }
               if(round_current_pos<=0)
               {
                open();
                   SqueezeB();
                   
                   Go_Counter();
                   gripB=1;
                   gripA=0;
                   return 1;
               }
           }
        }

        if(order_counter == order)
        {
            if(gripB==1)
           {
               if(round_current_pos>0)
               {
                   Go_Clock();
                   return 1;
               }
               if(round_current_pos<=0)
               {
                openB();
                   Squeeze();
                   
                   Go_Counter();
                   gripB=0;
                   gripA=1;
                   return 1;
               }
           }
           if(gripB==0)
           {
               if((max_rad-round_current_pos)>0)
               {
                   Go_Counter();
                   return 1;
               }
               if((max_rad-round_current_pos)<=0)
               {
                open();
                   SqueezeB();
                   
                   Go_Clock();
                   gripB=1;
                   gripA=0;
                   return 1;
               }
           }
        }
    }

    if(order!=round_last)//动作指令前后不一致
    {
        if(order_stop_round == order)//停止指令
        {
            return 1;
        }
        
        if (order_clock==order)
        {
            if(gripB==1)
           {
               if((max_rad-round_current_pos)>0)
               {
                   Go_Counter();
                   return 1;

               }
               if((max_rad-round_current_pos)<=0)
               {
                openB();
                   Squeeze();
                   
                   Go_back();
                   gripB=0;
                   gripA=1;
                   return 1;
               }
           }
           if(gripB==0)
           {
               if((max_rad-round_current_pos)>0)
               {
                open();
                   SqueezeB();
                   
                   gripB=1;
                   gripA=0;
                   Go_Counter();
                   return 1;
               }
               if((max_rad-round_current_pos)<=0)
               {
                   Go_Clock();
                   return 1;
               }
           }
        }
        if (order_counter==order)
        {
            if(gripB==1)
           {
               if(round_current_pos>0)
               {
                   Go_Clock();
                   return 1;
               }
               if(round_current_pos<=0)
               {
                openB();
                  Squeeze();
                   
                   Go_back();
                   gripB=0;
                   gripA=1;
                   return 1;
               }
           }
           if(gripB==0)
           {
               if(round_current_pos>0)
               {
                open();
                   SqueezeB();
                   
                   gripB=1;
                   gripA=1;
                   return 1;
               }
               if(round_current_pos<=0)
               {
                   Go_Counter();
                   return 1;
               }
           }
        }
    }
    
}

void read_package()
{
    if(Line_Order==1)
      {
          Line_Order = order_back;
      }
      else if(Line_Order==2)
      {
          Line_Order = order_go;
      }
      else
      {
          Line_Order = order_stop_line;
      }




    if(Round_Order==1)
      {
          Round_Order = order_clock;
      }
      else if(Round_Order==2)
      {
          Round_Order = order_counter;
      }
      else
      {
          Round_Order = order_stop_round;
      }

}

void Go_ahead()//脉冲控制电机
{
    digitalWrite(dirPin,LOW);
    for(int i=0 ;i<STEPS_PER_REV;i++)
    {
        digitalWrite(stepPin,HIGH);
        delayMicroseconds(X_speed);
        digitalWrite(stepPin,LOW);
        delayMicroseconds(X_speed);
        line_current_pos++;
    }
}

void Go_back()
{
    digitalWrite(dirPin,HIGH);
    for(int i=0 ;i<STEPS_PER_REV;i++)
    {
        digitalWrite(stepPin,HIGH);
        delayMicroseconds(X_speed);
        digitalWrite(stepPin,LOW);
        delayMicroseconds(X_speed);
        line_current_pos--;
    }
}

void Go_Clock()
{
    digitalWrite(dirPin_2,HIGH);
    for(int i=0 ;i<STEPS_PER_REV;i++)
    {
        digitalWrite(stepPin_2,HIGH);
        delayMicroseconds(Y_speed);
        digitalWrite(stepPin_2,LOW);
        delayMicroseconds(Y_speed);
        round_current_pos--;
    }
}

void Go_Counter()
{
    digitalWrite(dirPin_2,LOW);
    for(int i=0 ;i<STEPS_PER_REV;i++)
    {
        digitalWrite(stepPin_2,HIGH);
        delayMicroseconds(Y_speed);
        digitalWrite(stepPin_2,LOW);
        delayMicroseconds(Y_speed);
        round_current_pos++;
    }
}


void Squeeze()//夹紧
{
    tiaozheng();
}

void open()//夹抓放松
{

    digitalWrite(A_chong,LOW);
    digitalWrite(A_pai,HIGH);
    delay(3000);//延时3s把气放完
    digitalWrite(A_pai,LOW);

}


void SqueezeB()//夹紧
{
    digitalWrite(B_pai,LOW);
    digitalWrite(B_chong,HIGH);
    delay(2500);//持续2.5s充气
    digitalWrite(B_chong,LOW);
    delay(500);
}

void openB()//夹抓放松
{
    digitalWrite(B_chong,LOW);
    digitalWrite(B_pai,HIGH);
    delay(3000);
    digitalWrite(B_pai,LOW);

}











void Sachong()
{
    digitalWrite(A_pai,LOW);
    digitalWrite(A_chong,HIGH);
    delay(25);
    digitalWrite(A_chong,LOW);
    delay(700);
}

void Sapai()
{
    digitalWrite(A_chong,LOW);
    digitalWrite(A_pai,HIGH);
    delay(25);
    digitalWrite(A_pai,LOW);
    delay(100);
}

void Sbchong()
{
    digitalWrite(B_pai,LOW);
    digitalWrite(B_chong,HIGH);
    delay(25);
    digitalWrite(B_chong,LOW);
    delay(100);
}

void Sbpai()
{
    digitalWrite(B_chong,LOW);
    digitalWrite(B_pai,HIGH);
    delay(25);
    digitalWrite(B_pai,LOW);
    delay(100);
}
