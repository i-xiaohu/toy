import socket


def start_tcp_client(ip, port):
    # create socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    while True:
        try:
            print("pstat-client: start connect to server")
            s.connect((ip, port))
            break
        except socket.error:
            print("Fail to connect to server")
            return

    # send and receive
    while True:
        print("connect success")
        count = 0
        while True:
            msg = 'This is %d times sending' % count
            s.send(msg.encode('utf-8'))
            print("send len is : [%d]" % len(msg))

            msg = s.recv(1024)
            print(msg.decode('utf-8'))
            print("recv len is : [%d]" % len(msg))

            count += 1

            if count == 5:
                msg = 'server_finish'
                print("total send times is : %d " % count)
                s.send(msg.encode('utf-8'))
                break
        break

    s.close()


if __name__ == '__main__':
    ip = socket.gethostbyname("login1")
    print(ip)
    ip = "192.168.1.23"
    start_tcp_client(ip, 6000)
