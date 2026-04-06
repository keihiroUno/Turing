#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define N 50
#define M 50
#define Nx 256
#define My 256
#define Tc 1e5

int main(void){
    FILE *data,*gp;
    double u[Nx+1][My+1],v[Nx+1][My+1];
    double u_new[Nx+1][My+1],v_new[Nx+1][My+1];
    double delta_x,delta_y,delta_t,alpha_x,alpha_y;
    double a=0.2,b=1.46;
    double beta=0;
    double us,vs;
    double difMax,max;
    double eps=1e-10;
    double d=20;
    int i,j,n,iu,id,ju,jd;
    int k;
    int timeskip=Tc;
    char result[256];

    us=(a+1)/b;
    vs=((a+1)/b)*((a+1)/b);

    delta_x=N/(double)Nx;
    delta_y=M/(double)My;
    delta_t=1e-5;
    alpha_x=delta_t/(delta_x*delta_x);
    alpha_y=delta_t/(delta_y*delta_y);
    printf("d*alpha_x+d*alpha_y=%f\n",d*alpha_x+d*alpha_y);

    n=0;

    //初期値の設定
    srand((unsigned)time(NULL));
    for(i=0;i<=Nx;i++){
        for(j=0;j<=My;j++){
            u[i][j]=us+((rand()*20.0/(1.0+RAND_MAX))-5.0)/100;
            v[i][j]=vs+((rand()*20.0/(1.0+RAND_MAX))-5.0)/100;
        }
    }

    //結果の抽出とファイル出力
    sprintf(result,"./Turing2_a02b146/data%d.txt",n);
    data=fopen(result,"w");
    if(data==NULL){
        printf("ファイルを開くことができませんでした: %s\n",result);
        return 1;
    }
    for(i=0;i<=Nx;i++){
        for(j=0;j<=My;j++){
            fprintf(data,"%f %f %.10f %.10f\n",i*delta_x,j*delta_y,u[i][j],v[i][j]);
        }
        fprintf(data,"\n");
    }
    fclose(data);

    //差分式の計算
    do{
        n++;
        for(i=0;i<=Nx;i++){
            iu=i+1;
            if(iu==Nx+1) iu=Nx-1;
            id=i-1;
            if(id==-1) id=1;
            for(j=0;j<=My;j++){
                ju=j+1;
                if(ju==My+1) ju=My-1;
                jd=j-1;
                if(jd==-1) jd=1;
                u_new[i][j] = u[i][j] + alpha_x*(u[iu][j]-2.0*u[i][j]+u[id][j]) + alpha_y*(u[i][ju]-2.0*u[i][j]+u[i][jd]) + delta_t*(a - b*u[i][j] + u[i][j]*u[i][j]/v[i][j]);
                v_new[i][j] = v[i][j] + d*alpha_x*(v[iu][j]-2.0*v[i][j]+v[id][j]) + d*alpha_y*(v[i][ju]-2.0*v[i][j]+v[i][jd]) + delta_t*(u[i][j]*u[i][j] - v[i][j] + beta);
            }
        }
        //収束判定
        difMax=0.0;
        for(i=0;i<=Nx;i++){
            for(j=0;j<=My;j++){
                max=(fabs(u_new[i][j]-u[i][j])+fabs(v_new[i][j]-v[i][j]))/2.0;
                if(max>difMax) difMax=max;
            }
        }

        if(n%timeskip==0){
            sprintf(result,"./Turing2_a02b146/data%d.txt",n);
            data=fopen(result,"w");
            for(i=0;i<=Nx;i++){
                for(j=0;j<=My;j++){
                    fprintf(data,"%f %f %.10f %.10f\n",i*delta_x,j*delta_y,u_new[i][j],v_new[i][j]);
                }
                fprintf(data,"\n");
            }
            fclose(data);
            printf("%d %.15f\n",n,difMax);
        }
        for(j=0;j<=My;j++){
            for(i=0;i<=Nx;i++){
                u[i][j]=u_new[i][j];
                v[i][j]=v_new[i][j];
            }
        }

        // 途中で止める
        if(n==1e8) break;
            
    }while(difMax>eps || n%timeskip!=0);
    printf("%f\n",a);
    return 0;
}
