using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.IO;
using System.Threading;
using System.Diagnostics;

public class button : MonoBehaviour
{
    void Start()
    {
        Button btn = this.GetComponent<Button>();
        btn.onClick.AddListener(OnClick);
    }

    private void OnClick()
    {
        UnityEngine.Debug.Log("Button Clicked. ClickHandler.");
        float[] va = new float[26];
        string[] vastr = new string[26];
        for(int i = 0; i < 26; ++ i)
        {
            string name = "Slider " + (i + 1);
            va[i] = GameObject.Find(name).GetComponent<Slider>().value;
            vastr[i] = va[i].ToString();
        }
        string filename = Application.dataPath + "/Resources/measure.txt";
        File.WriteAllLines(filename, vastr);

        Thread newThread = new Thread(new ThreadStart(RunShellThreadStart));
        newThread.Start();

    }


    private static void RunShellThreadStart()
    {
        string cmdTxt = "python C:\\Users\\amomorning\\source\\repos\\semantic_body\\bodyViz\\Assets\\Codes\\predict.py";

        RunCommand(cmdTxt);

        string cmdtt = "C:\\Users\\amomorning\\source\\repos\\semantic_body\\x64\\Release\\postprocessing.exe";
        RunCommand(cmdtt);
        //RunProcessCommand("notepad", @"C:\Users\pc\Desktop\1.txt");
        //RunProcessCommand("explorer.exe", @"D:\");
        UnityEngine.Debug.Log("Success");

    }

    private static void RunCommand(string command)
    {
        Process process = new Process();
        process.StartInfo.FileName = "powershell";
        process.StartInfo.Arguments = command;

        process.StartInfo.CreateNoWindow = false; // 获取或设置指示是否在新窗口中启动该进程的值（不想弹出powershell窗口看执行过程的话，就=true）
        process.StartInfo.ErrorDialog = false; // 该值指示不能启动进程时是否向用户显示错误对话框
        process.StartInfo.UseShellExecute = false;
        //process.StartInfo.RedirectStandardError = true;
        //process.StartInfo.RedirectStandardInput = true;
        //process.StartInfo.RedirectStandardOutput = true;

        process.Start();

        //process.StandardInput.WriteLine(@"explorer.exe D:\");
        //process.StandardInput.WriteLine("pause");

        process.WaitForExit();
        process.Close();
    }


}
